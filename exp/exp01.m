
exp_id = 'exp01';
exp_description = 'This is the first experiment ever.';

machines = {'bajel3am@u-pl28', 'bajel3am@u-pl29'};

% BBOB parameters
bbParamDef(1).name   = 'dimensions';
bbParamDef(1).values = {2, 5 10};
bbParamDef(2).name   = 'functions';
bbParamDef(2).values = {2, 3, 8};
% dimensions  = [10];     % which dimensions to optimize, subset of [2 3 5 10 20 40];
% functions   = [8];      % function ID's to optimize (2 Sphere, 3 Rastrigin, 8 Rosenbrock)
bbParamDef(3).name   = 'opt_function';
bbParamDef(3).values = {@opt_s_cmaes};
% opt_function = @opt_s_cmaes;    % function being optimized -- BBOB wrap-around with header
%                                 % xbest = function( fun, dim, ftarget, maxfunevals )

% Surrogate model parameter lists
sgParamDef(1).name   = 'evoControl';
sgParamDef(1).values = {'individual'};
sgParamDef(2).name   = 'modelType';
sgParamDef(2).values = {'gp', 'rf'};
sgParamDef(3).name   = 'evoControlPreSampleSize';
sgParamDef(3).values = {0.25, 0.5, 0.75};
sgParamDef(4).name   = 'evoControlIndividualExtension';
sgParamDef(4).values = {10};
sgParamDef(5).name   = 'evoControlBestFromExtension';
sgParamDef(5).values = {0.1};
sgParamDef(6).name   = 'evoControlTrainRange';
sgParamDef(6).values = {2, 3, 4, 6, 8};

% CMA-ES parameters
cmParamDef(1).name   = 'PopSize';
cmParamDef(1).values = {'(4 + floor(3*log(N)))', '(8 + floor(6*log(N)))'};


% Divide instances to machines
params = [bbParams, sgParams, cmParams];
nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
nMachines = length(machines);

combsPerMachine = ceil(nCombinations / nMachines);
startIdxs = 1:combsPerMachine:nCombinations;
endIdxs =   0:combsPerMachine:(nCombinations-1);




% BBOB constant parameters
instances = [1:5, 31:40];       % [1:5, 31:40]
maxfunevals = '250 * dim';      % 10*dim is a short test-experiment taking a few minutes 
                                % INCREMENT maxfunevals successively to larger value(s)
minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm
bbobpath = 'vendor/bbob';    % should point to fgeneric.m etc.
datapath = ['../log/bbob/' exp_id];  % different folder for each experiment
  pathstr = fileparts(mfilename('fullpath'));
  datapath = [pathstr filesep datapath];
  addpath([pathstr filesep bbobpath]);
opt.algName = exp_description;
opt.comments = '';

% runs an experiment for benchmarking MY_OPTIMIZER
% on the noise-free testbed. fgeneric.m and benchmarks.m
% must be in the path of Matlab/Octave

more off;  % in octave pagination is on by default

t0 = clock;
rand('state', sum(100 * t0));

results = cell(0);
y_evals = cell(0);

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

    for iinstance = instances   % 15 function instances
      fgeneric('initialize', ifun, iinstance, datapath, opt); 

      % independent restarts until maxfunevals or ftarget is reached
      for restarts = 0:maxrestarts
        if restarts > 0  % write additional restarted info
          fgeneric('restart', 'independent restart')
        end
        [xopt, ilaunch, ye] = bbParams.opt_function('fgeneric', dim, fgeneric('ftarget'), ...
                     eval(maxfunevals) - fgeneric('evaluations'));
        % we don't have this information from CMA-ES :(
        % results = cat(1,results,res);
        % ye = [res.deltasY res.evaluations];
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

      fgeneric('finalize');
      save([datapath filesep exp_id '_tmp.mat'], 'results', 'y_evals');
    end
    disp(['      date and time: ' num2str(clock, ' %.0f')]);

    elapsedTime = toc(t);

    exp_results.evals = inst_results_evals;
    exp_results.restarts = inst_results_restarts;
    exp_results.y_evals = y_evals;
    exp_results.time = elapsedTime;

    save([datapath filesep exp_id '_' num2str(ifun) '_' num2str(dim) 'D.mat'], 'exp_id', 'exp_settings', 'exp_results', 'results', 'y_evals');

  end
  disp(sprintf('---- dimension %d-D done ----', dim));
end


