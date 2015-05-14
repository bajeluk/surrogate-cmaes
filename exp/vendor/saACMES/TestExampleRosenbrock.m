global settings;

settings.instances = [1];   % instance of the problem
settings.dims = [10];       % problem dimension
settings.funs = [8];        % problem/function id in BBOB framework
settings.pathname = 'd1';
settings.algname = '';
settings.ntarray = [1];
settings.savfile = 'r1';

settings.BIPOP = 0;         % BIPOP if BIPOP = 1 and IPOP otherwise
settings.newRestartRules = 0;
settings.noisy = 0;         % chagne few parameters if problem is noisy (noisy = 1)
settings.CMAactive = 1;     % use active CMA or not, active usually performs better
settings.withFileDisp = 1;  % save some information if files?
settings.withSurr = 1;      % use surrogate models? if withSurr=0, than we have original IPOP or BIPOP-CMA-ES
settings.modelType = 1;
settings.withModelEnsembles = 0;
settings.withModelOptimization = 0; % optimize model hyper-parameters or we know some good default parameters? ( see xacmes.m ) =0 is cheap, =1 is 20 times (here) more expensive
settings.hyper_lambda = 20; % if optimize hyper-parameters, how much models build each iteration? (population size of the second CMA-ES)
settings.iSTEPminForHyperOptimization = 1;  % if model control says that model should be used for less than iSTEPminForHyperOptimization generation, 
                                            % than we will not even try to optimize hyper-parameters and spend the computational effort on that
settings.MaxEvals = '1e6*dim';  % maximum budget of function evaluation can be used during the optimization
settings.MaxEvalsWithSurrogate = '1e4*20';
settings.lambdaMult = 1;
settings.muMult = 1;
settings.largeLambdaMinIter = 3;


settings.withDisp = 0;          % show some graphs and intermediate results?
settings.maxStepts = 20;        % maximum number of generation the surrogate model can be used instead of real function
settings.maxerr = 0.45;         % see details in the paper
settings.alpha = 0.20;          % see details in the paper
settings.iterstart = 10;        % starting from iterstart iteration we will use surrogate models

Adapter();                      % launch the algorithm on BBOB problem
                          
% to launch the algorithm on your own function, change MyFunc.m and all feval('fgeneric',..), fgeneric(...) you can find in the code