function bbob_test_01(id, exp_id, path)
  exppath = [path filesep exp_id];
  load([exppath filesep 'scmaes_params.mat']);
  [bbParams, surrogateParams, cmaesParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  addpath(exppath);
  
  % BBOB constant parameters
  minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
  maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm
  bbobpath = 'vendor/bbob';    % should point to fgeneric.m etc.
  datapath = ['../log/bbob/' exp_id];  % different folder for each experiment
    pathstr = fileparts(mfilename('fullpath'));
    datapath = [pathstr filesep datapath];
    addpath([pathstr filesep bbobpath]);
  % opt.algName = exp_description;
  opt.comments = '';

  % runs an experiment for benchmarking MY_OPTIMIZER
  % on the noise-free testbed. fgeneric.m and benchmarks.m
  % must be in the path of Matlab/Octave

  more off;  % in octave pagination is on by default

  t0 = clock;
  rand('state', sum(100 * t0));

  results = cell(0);
  y_evals = cell(0);
  instances = bbParams.instances;
  maxfunevals = bbParams.maxfunevals;

  for dim = bbParams.dimensions            % small dimensions first, for CPU reasons
    % for ifun = benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
    for ifun = bbParams.functions          % or benchmarksnoisy(...)

      t = tic;
      exp_settings.dim = dim;
      exp_settings.bbob_function = ifun;
      exp_settings.exp_id = exp_id;
      exp_settings.instances = instances;
      inst_results_evals = [];
      inst_results_restarts = [];
      inst_results_fbests = [];
      inst_results_f025 = [];
      inst_results_f050 = [];
      inst_results_f075 = [];
      inst_results_stopflags = {};
      fmin = Inf;
      evalsRestartCorrection = 0;

      for iinstance = instances   % 15 function instances
        fgeneric('initialize', ifun, iinstance, datapath, opt); 

        % independent restarts until maxfunevals or ftarget is reached
        for restarts = 0:maxrestarts
          if restarts > 0  % write additional restarted info
            fgeneric('restart', 'independent restart')
          end
          [xopt, ilaunch, ye, stopflag] = bbParams.opt_function('fgeneric', dim, fgeneric('ftarget'), ...
                      eval(maxfunevals) - fgeneric('evaluations'), id);
          % we don't have this information from CMA-ES :(
          % results = cat(1,results,res);
          % ye = [res.deltasY res.evaluations];

          if (fmin < Inf)
            ye(:,1) = min([ye(:,1) repmat(fmin,size(ye,1),1)], [], 2);
            ye(:,2) = ye(:,2) + evalsRestartCorrection;
          end
          fmin = min([ye(:,1); fmin]);
          evalsRestartCorrection = fgeneric('evaluations');
          y_evals = cat(1,y_evals,ye);

          if fgeneric('fbest') < fgeneric('ftarget') || ...
            fgeneric('evaluations') + eval(minfunevals) > eval(maxfunevals)
            break;
          end  
        end

        disp(sprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                      ' fbest-ftarget=%.4e, elapsed time [h]: %.2f'], ...
                    ifun, dim, iinstance, ...
                    fgeneric('evaluations'), ...
                    restarts, ...
                    fgeneric('fbest') - fgeneric('ftarget'), ...
                    etime(clock, t0)/60/60));

        inst_results_evals = [inst_results_evals fgeneric('evaluations')];
        inst_results_restarts = [inst_results_restarts restarts];
        inst_results_fbests = [inst_results_fbests min(y_evals{end}(:,1))];
        inst_results_f025   = [inst_results_f025 y_evals{end}( floor(size(y_evals{end},1)/4) ,1)];
        inst_results_f050   = [inst_results_f050 y_evals{end}( floor(size(y_evals{end},1)/2) ,1)];
        inst_results_f075   = [inst_results_f075 y_evals{end}( floor(3*size(y_evals{end},1)/4) ,1)];
        inst_results_stopflags{end+1} = stopflag;

        fgeneric('finalize');
        tmpFile = [exppath filesep exp_id '_' num2str(id) '_tmp.mat'];
        save(tmpFile, 'results', 'y_evals');
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

      save([exppath filesep exp_id '_' num2str(ifun) '_' num2str(dim) 'D_' num2str(id) '.mat'], 'exp_id', 'exp_settings', 'exp_results', 'results', 'y_evals', 'surrogateParams', 'cmaesParams');

      delete(tmpFile);
    end
    disp(sprintf('---- dimension %d-D done ----', dim));
  end
end
