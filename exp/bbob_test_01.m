function bbob_test_01(id, exp_id, exppath_short, varargin)
%BBOB_TEST_01 Run BBOB experiment with parameters defined by number 'id'
% parameters:
%   id          the number which correspond to one set of parameters
%               (BBOB parameters, surrogateOpts and CMA-ES parameters)
%   exp_id      unique string identifier of the experiment
%   exppath_short        directory where experiment output data will be placed
%   OPTARG1     directory where local files should be placed ("scratchdir" at metacentrum)

  % GNUPlot script where special strings will be replaced
  gnuplotScript = 'twoAlgsPlotExtended.gpi';

  % GnuPlot script should be in $ALGROOT/exp/
  gnuplotScript = [exppath_short filesep '..' filesep gnuplotScript];
  % Directory for internal results of _this_ function
  exppath = [exppath_short filesep exp_id];
  load([exppath filesep 'scmaes_params.mat']);
  [bbParams, surrogateParams, cmaesParams, nNonBbobValues] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  
  % BBOB parameters
  minfunevals = 'dim + 2';      % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
  bbobpath = 'vendor/bbob';     % should point to fgeneric.m etc.
  maxrestarts = defopts(bbParams, 'maxrestarts', 1e4); % SET to zero for an entirely deterministic algorithm 
  % PROGRESS_LOG = defopts(bbParams, 'progressLog', false); % surrogateModel progress log

  localDatapath = [];       % directory in the shared folder where results of each instance will be copied through the progress
  if (nargin >= 4 && ~isempty(varargin{1}))
    % datapath is typically on the computing node in   $SCRATCHDIR/job_output/bbob_output:
    datapath = [varargin{1} filesep 'bbob_output'];
    if (isempty(strfind(datapath, exppath_short)))
      % localDatapath is typically on the NFS server in
      %   $HOME/prg/surrogate-cmaes/exp/experiments/$EXPID/bbob_output_tmp
      % and the BBOB results are stored here after each completed instance
      localDatapath = [exppath filesep 'bbob_output_tmp'];
      [~, ~] = mkdir(localDatapath);
    end
  else
    datapath = [exppath filesep 'bbob_output'];
  end
  [~, ~] = mkdir(datapath);

  % opt.algName = exp_description;
  opt.comments = '';

  more off;  % in octave pagination is on by default

  t0 = clock;
  % Initialize random number generator
  exp_settings.seed = myeval(defopts(bbParams, 'seed', 'floor(sum(100 * clock()))'));
  rng(exp_settings.seed);

  instances = bbParams.instances;
  maxfunevals = bbParams.maxfunevals;

  try

  for dim = bbParams.dimensions            % small dimensions first, for CPU reasons
    % for ifun = benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
    for ifun = bbParams.functions          % or benchmarksnoisy(...)

      % ===== OUR INTERESTING RESULTS =====

      exp_settings.dim = dim;
      exp_settings.bbob_function = ifun;
      exp_settings.exp_id = exp_id;
      exp_settings.instances = instances;
      exp_settings.resume = defopts(bbParams, 'resume', false);
      exp_settings.progressLog = defopts(bbParams, 'progressLog', false);

      expFileID = [num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
      resultsFile = [exppath filesep exp_id '_results_' expFileID];
      opt.algName = [exp_id '_' expFileID];
      datapath = [datapath filesep expFileID];
      [~, ~] = mkdir(datapath);
      cmaes_out = [];

      [exp_results, tmpFile, cmaes_out] = runTestsForAllInstances(bbParams.opt_function, id, exp_settings, datapath, opt, maxrestarts, eval(maxfunevals), eval(minfunevals), t0, exppath, localDatapath, false);

      y_evals = exp_results.y_evals;

      save([resultsFile '.mat'], 'exp_id', 'exp_settings', 'exp_results', 'y_evals', 'surrogateParams', 'cmaesParams', 'bbParams', 'cmaes_out');

      % ===== PURE CMAES RESULTS =====

      % Tripple the number of pure CMA-ES instances
      exp_settings.instances = [instances instances+50 instances+100];

      cmaesId = floor((id-1) / nNonBbobValues) * nNonBbobValues + 1;
      % test if pure CMA-ES results exist; if no, generate them
      cmaesResultsFile = [exppath filesep 'cmaes_results' filesep exp_id '_purecmaes_' num2str(ifun) '_' num2str(dim) 'D_' num2str(cmaesId) '.mat'];
      if (~ exist(cmaesResultsFile, 'file'))
        opt.algName = [exp_id '_' expFileID '_cmaes'];
        exp_cmaes_results = runTestsForAllInstances(@opt_cmaes, id, exp_settings, datapath, opt, maxrestarts, eval(maxfunevals), eval(minfunevals), t0, exppath, localDatapath, true);

        % test if the results still doesn't exist, if no, save them :)
        if (~ exist(cmaesResultsFile, 'file'))
          y_evals = exp_cmaes_results.y_evals;
          save(cmaesResultsFile, 'exp_id', 'exp_settings', 'exp_cmaes_results', 'y_evals', 'surrogateParams', 'cmaesParams');
        end
      end

      % Set the number of instances back to the original value
      exp_settings.instances = instances;

      % ===== GENERATE PICTURE =====

      % Load the pure CMA-ES results
      load(cmaesResultsFile, 'exp_cmaes_results');

      % Save the data for gnuplot
      gnuplotFile = [exppath filesep exp_id '_gnuplot_' num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
      generateGnuplotDataExtended([gnuplotFile '.dat'], exp_results, exp_cmaes_results, eval(maxfunevals));

      % save gnuplot script
      if (isfield(surrogateParams, 'modelType')); modelType = surrogateParams.modelType;
      else modelType = ''; end
      if (isfield(surrogateParams, 'evoControl')); evoControl = surrogateParams.evoControl;
      else evoControl = ''; end

      gnuplotScriptCommand = ['sed "s#\<DATAFILE\>#' gnuplotFile '.dat#; s#\<OUTPUTFILE\>#' resultsFile '#; s#\<TITLE\>#f' num2str(ifun) ', ' num2str(dim) 'D#; s#\<DATALINETITLE\>#' modelType ' surrogate, ' evoControl ' EC#; s#\<PARAMS1\>#' sprintfStruct(surrogateParams, 'escape') '#; s#\<PARAMS2\>#' sprintfStruct(exp_settings, 'escape') '#" ' gnuplotScript ' > ' gnuplotFile '.gpi'];
      disp(gnuplotScriptCommand);
      system(gnuplotScriptCommand);
      % call gnuplot
      system(['LD_LIBRARY_PATH="" gnuplot ' gnuplotFile '.gpi']);

      % print out settings into the text-file
      fid = fopen([resultsFile '.txt'], 'w');
      printSettings(fid, exp_settings, exp_results, surrogateParams, cmaesParams);
      fclose(fid);

      delete(tmpFile);
    end
    fprintf('---- dimension %d-D done ----\n', dim);
  end

  catch err
    save([resultsFile '_ERROR.mat']);
    fprintf('#########################################################\n');
    fprintf('#########################################################\n');
    fprintf('              Matlab ended with error!\n');
    fprintf('---------------------------------------------------------\n');
    fprintf('%s\n', err.identifier);
    fprintf('%s\n', err.message);
    fprintf('---------------------------------------------------------\n');
    getReport(err)
    fprintf('---------------------------------------------------------\n');
    if (exist('exp_results', 'var'))
      fprintf('---------------------------------------------------------\n');
      printSettings(1,  exp_settings, exp_results, surrogateParams, cmaesParams);
    end
    fprintf('#########################################################\n');
    fprintf('#########################################################\n');
    % comment the following "exit(1)" when debugging -- it shutdowns the
    % whole Matlab if an error occures
    exit(1);
    throw(err);
  end
end

function [exp_results, tmpFile, cmaes_out] = runTestsForAllInstances(opt_function, id, exp_settings, datapath, opt, maxrestarts, maxfunevals, minfunevals, t0, exppath, localDatapath, isPureCmaes)
  y_evals = cell(0);
  cmaes_out = cell(0);

  exp_results.evals = [];
  exp_results.restarts = [];
  exp_results.fbests = [];
  exp_results.f025 = [];
  exp_results.f050 = [];
  exp_results.f075 = [];
  exp_results.stopflags = {};
  exp_results.time = 0;
  evalsRestartCorrection = 0;

  tmpFile = [exppath filesep exp_settings.exp_id '_tmp_' num2str(id) '.mat'];

  % load interrupted "_tmp" results if exp_settings.resume is set
  [datapathRoot, expFileID] = fileparts(datapath);
  if (~isPureCmaes && exp_settings.resume && ~isempty(localDatapath) ...
      && exist([localDatapath filesep expFileID], 'dir') ...
      && exist(tmpFile, 'file'))
    [nCompletedInstances, y_evals, exp_results, cmaes_out] = loadInterruptedInstances(tmpFile);
    system(['cp -pR ' localDatapath '/' expFileID ' ' datapathRoot]);
    % copy also not-finished logs of this experiment ID
    % TODO: test this!
    system(['cp -pR ' localDatapath '/' exp_settings.exp_id '_log_' expFileID '.dat ' datapathRoot]);
  else
    nCompletedInstances = 0;
  end

  for iinstance = exp_settings.instances((nCompletedInstances+1):end)   % 15 function instances
    fmin = Inf;

    fgeneric('initialize', exp_settings.bbob_function, iinstance, datapath, opt); 
    yeRestarts = [];
    cmaes_out{end+1}  = {};
    t = tic;
    if (~isPureCmaes && isfield(exp_results, 'rngState'));
      fprintf('Loading saved random number generator state (%d, %d)...\n', exp_results.rngState.Seed, exp_results.rngState.State(1));
      rng(exp_results.rngState);
    end
    xstart = 8 * rand(exp_settings.dim, 1) - 4;
    restartMaxfunevals = maxfunevals;

    % independent restarts until maxfunevals or ftarget is reached
    for restarts = 0:maxrestarts
      if restarts > 0  % write additional restarted info
        fgeneric('restart', 'independent restart')
      end
      [xopt, ye, stopflag, cmaes_out_1] = opt_function('fgeneric', exp_settings.dim, fgeneric('ftarget'), ...
                  restartMaxfunevals, id, exppath, xstart, fileparts(datapath), iinstance);
      if (exp_settings.progressLog)
        cmaes_out{end}{end+1} = cmaes_out_1;
      end

      if (fmin < Inf)
        ye(:,1) = min([ye(:,1) repmat(fmin,size(ye,1),1)], [], 2);
        ye(:,2) = ye(:,2) + evalsRestartCorrection;
      end
      fmin = min([ye(:,1); fmin]);
      yeRestarts = [yeRestarts; ye];
      evalsRestartCorrection = fgeneric('evaluations');

      if fgeneric('fbest') < fgeneric('ftarget') || ...
        fgeneric('evaluations') + minfunevals > restartMaxfunevals
        break;
      else
        % try to improve the best foud solution
        xstart = xopt;
        restartMaxfunevals = restartMaxfunevals - fgeneric('evaluations');
      end  
    end

    elapsedTime = toc(t);
    y_evals = cat(1,y_evals,yeRestarts);

    fprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                  ' fbest-ftarget=%.4e, elapsed time [h]: %.2f\n'], ...
                exp_settings.bbob_function, exp_settings.dim, iinstance, ...
                fgeneric('evaluations'), ...
                restarts, ...
                fgeneric('fbest') - fgeneric('ftarget'), ...
                etime(clock, t0)/60/60);

    exp_results.evals(end+1)  = fgeneric('evaluations');
    exp_results.restarts(end+1) = restarts;
    exp_results.fbests(end+1) = min(y_evals{end}(:,1));
    exp_results.f025(end+1)   = y_evals{end}( max([1 floor(size(y_evals{end},1)/4)]) ,1);
    exp_results.f050(end+1)   = y_evals{end}( max([1 floor(size(y_evals{end},1)/2)]) ,1);
    exp_results.f075(end+1)   = y_evals{end}( max([1 floor(3*size(y_evals{end},1)/4)]) ,1);
    exp_results.stopflags{end+1} = stopflag;
    exp_results.y_evals       = y_evals;
    exp_results.time          = exp_results.time + elapsedTime;

    fgeneric('finalize');
    exp_id = exp_settings.exp_id;
    if (~isPureCmaes)
      exp_results.rngState = rng();
      save(tmpFile, 'exp_settings', 'exp_id', 'y_evals', 'exp_results', 'cmaes_out');
    end

    % copy the output to the final storage (if OUTPUTDIR and EXPPATH differs)
    if (~isempty(localDatapath) && isunix ...
        && (~isPureCmaes || iinstance == exp_settings.instances(end)) )
      system(['cp -pR ' datapath ' ' localDatapath '/']);
      system(['cp -pR ' datapath '/../' exp_settings.exp_id '_log_' expFileID '.dat ' localDatapath]);
    end
  end
  disp(['      date and time: ' num2str(clock, ' %.0f')]);
end

function printSettings(fid, exp_settings, exp_results, surrogateParams, cmaesParams)
  fprintf(fid, '===== Experiment: %s =====\n\n', exp_settings.exp_id);
  fprintf(fid, '== BBOB experiment settings: ==\n');
  fprintf(fid, sprintfStruct(exp_settings));
  fprintf(fid, '%15s: %f\n', 'time elapsed', exp_results.time);
  fprintf(fid, '\n== Surrogate model parameters: ==\n');
  fprintf(fid, sprintfStruct(surrogateParams));
  fprintf(fid, '\n== CMA-ES parameters: ==\n');
  fprintf(fid, sprintfStruct(cmaesParams));
  fprintf(fid, '\n== CMA-ES surrogate model options: ==\n');
  if (isfield(surrogateParams, 'modelOpts'))
    fprintf(fid, sprintfStruct(surrogateParams.modelOpts)); end
  fprintf(fid, '\n== Numerical results: ==\n\n');
  fprintf(fid, 'fbests:\n%s\n\n', num2str(exp_results.fbests));
  fprintf(fid, 'f075:\n%s\n\n', num2str(exp_results.f075));
  fprintf(fid, 'f050:\n%s\n\n', num2str(exp_results.f050));
  fprintf(fid, 'f025:\n%s\n\n', num2str(exp_results.f025));
end

function [nCompletedInstances, y_evals, exp_results, cmaes_out] = loadInterruptedInstances(tmpFile)
  fprintf('Resuming previously interrupted experiment run...\n');
  fprintf('Loading results from:  %s\n', tmpFile);
  load(tmpFile);
  nCompletedInstances = size(y_evals, 1);
  fprintf('Completed instances (%d):  %s\n', nCompletedInstances, num2str(exp_settings.instances(1:nCompletedInstances)));
end
