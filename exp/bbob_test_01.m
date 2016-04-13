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

  % surrogateModel progress log
  PROGRESS_LOG = 0;

  % GnuPlot script should be in $ALGROOT/exp/
  gnuplotScript = [exppath_short filesep '..' filesep gnuplotScript];
  % Directory for internal results of _this_ function
  exppath = [exppath_short filesep exp_id];
  load([exppath filesep 'scmaes_params.mat']);
  [bbParams, surrogateParams, cmaesParams, nNonBbobValues] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  
  % BBOB constant parameters
  minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
  maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm
  bbobpath = 'vendor/bbob';    % should point to fgeneric.m etc.
  % addpath cannot be in deployed code!
  % addpath([exppath_short filesep '..' filesep bbobpath]);
  localDatapath = [];       % directory in the shared folder where results of each instance will be copied through the progress
  if (nargin >= 4 && ~isempty(varargin{1}))
    datapath = [varargin{1} filesep 'bbob_output'];
    if (isempty(strfind(datapath, exppath_short)))
      localDatapath = [exppath filesep 'bbob_output_tmp'];
      [~, ~] = mkdir(localDatapath);
    end
  else
    datapath = [exppath filesep 'bbob_output'];
    % old version:
    % datapath = ['../log/bbob/' exp_id];  % different folder for each experiment
  end
  [~, ~] = mkdir(datapath);

  % opt.algName = exp_description;
  opt.comments = '';

  % runs an experiment for benchmarking MY_OPTIMIZER
  % on the noise-free testbed. fgeneric.m and benchmarks.m
  % must be in the path of Matlab/Octave

  more off;  % in octave pagination is on by default

  t0 = clock;
  rand('state', sum(100 * t0));

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

      expFileID = [num2str(ifun) '_' num2str(dim) 'D_' num2str(id)];
      resultsFile = [exppath filesep exp_id '_results_' expFileID];
      opt.algName = [exp_id '_' expFileID];
      datapath = [datapath filesep expFileID];
      [~, ~] = mkdir(datapath);
      cmaes_out = [];

      if (PROGRESS_LOG)
        [exp_results, tmpFile, cmaes_out] = runTestsForAllInstances(bbParams.opt_function, id, exp_settings, datapath, opt, maxrestarts, eval(maxfunevals), eval(minfunevals), t0, exppath, localDatapath);
      else
        [exp_results, tmpFile] = runTestsForAllInstances(bbParams.opt_function, id, exp_settings, datapath, opt, maxrestarts, eval(maxfunevals), eval(minfunevals), t0, exppath, localDatapath);
      end

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
        exp_cmaes_results = runTestsForAllInstances(@opt_cmaes, id, exp_settings, datapath, opt, maxrestarts, eval(maxfunevals), eval(minfunevals), t0, exppath, '');

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
      system(['gnuplot ' gnuplotFile '.gpi']);

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
    getReport(err);
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

function [exp_results, tmpFile, cmaes_out] = runTestsForAllInstances(opt_function, id, exp_settings, datapath, opt, maxrestarts, maxfunevals, minfunevals, t0, exppath, localDatapath)
  y_evals = cell(0);
  cmaes_out = cell(0);

  t = tic;
  inst_results_evals = [];
  inst_results_restarts = [];
  inst_results_fbests = [];
  inst_results_f025 = [];
  inst_results_f050 = [];
  inst_results_f075 = [];
  inst_results_stopflags = {};
  evalsRestartCorrection = 0;

  for iinstance = exp_settings.instances   % 15 function instances
    fmin = Inf;

    fgeneric('initialize', exp_settings.bbob_function, iinstance, datapath, opt); 
    yeRestarts = [];

    % independent restarts until maxfunevals or ftarget is reached
    for restarts = 0:maxrestarts
      if restarts > 0  % write additional restarted info
        fgeneric('restart', 'independent restart')
      end
      [xopt, ilaunch, ye, stopflag, cmaes_out_1] = opt_function('fgeneric', exp_settings.dim, fgeneric('ftarget'), ...
                  maxfunevals, id, exppath);
      % we don't have this information from CMA-ES :(
      % ye = [res.deltasY res.evaluations];
      cmaes_out{end+1} = cmaes_out_1;

      if (fmin < Inf)
        ye(:,1) = min([ye(:,1) repmat(fmin,size(ye,1),1)], [], 2);
        ye(:,2) = ye(:,2) + evalsRestartCorrection;
      end
      fmin = min([ye(:,1); fmin]);
      yeRestarts = [yeRestarts; ye];
      evalsRestartCorrection = fgeneric('evaluations');

      if fgeneric('fbest') < fgeneric('ftarget') || ...
        fgeneric('evaluations') + minfunevals > maxfunevals
        break;
      end  
    end

    y_evals = cat(1,y_evals,yeRestarts);

    fprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                  ' fbest-ftarget=%.4e, elapsed time [h]: %.2f\n'], ...
                exp_settings.bbob_function, exp_settings.dim, iinstance, ...
                fgeneric('evaluations'), ...
                restarts, ...
                fgeneric('fbest') - fgeneric('ftarget'), ...
                etime(clock, t0)/60/60);

    inst_results_evals = [inst_results_evals fgeneric('evaluations')];
    inst_results_restarts = [inst_results_restarts restarts];
    inst_results_fbests = [inst_results_fbests min(y_evals{end}(:,1))];
    inst_results_f025   = [inst_results_f025 y_evals{end}( max([1 floor(size(y_evals{end},1)/4)]) ,1)];
    inst_results_f050   = [inst_results_f050 y_evals{end}( max([1 floor(size(y_evals{end},1)/2)]) ,1)];
    inst_results_f075   = [inst_results_f075 y_evals{end}( max([1 floor(3*size(y_evals{end},1)/4)]) ,1)];
    inst_results_stopflags{end+1} = stopflag;

    fgeneric('finalize');
    tmpFile = [exppath filesep exp_settings.exp_id '_tmp_' num2str(id) '.mat'];
    exp_id = exp_settings.exp_id;
    save(tmpFile, 'exp_settings', 'exp_id', 'y_evals');

    % copy the output to the final storage (if OUTPUTDIR and EXPPATH differs)
    if (~isempty(localDatapath) && isunix)
      system(['cp -pR ' datapath ' ' localDatapath '/']);
    end
  end
  disp(['      date and time: ' num2str(clock, ' %.0f')]);

  elapsedTime = toc(t);

  exp_results.evals = inst_results_evals;
  exp_results.restarts = inst_results_restarts;
  exp_results.f025 = inst_results_f025;
  exp_results.f050 = inst_results_f050;
  exp_results.f075 = inst_results_f075;
  exp_results.fbests = inst_results_fbests;
  exp_results.stopflags = inst_results_stopflags;
  exp_results.y_evals = y_evals;
  exp_results.time = elapsedTime;
end
