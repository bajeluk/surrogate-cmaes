global settings;

settings.instances = [5];
settings.dims = [20];
settings.funs = [23,24,3,4];
settings.pathname = 'd50';
settings.algname = '';
settings.ntarray = [1];
settings.savfile = 'r50';

settings.BIPOP = 1; 
settings.newRestartRules = 0; 
settings.noisy = 0;
settings.CMAactive = 1;
settings.withFileDisp = 1;
settings.withSurr = 1;
settings.modelType = 1;
settings.withModelEnsembles = 0;
settings.withModelOptimization = 1;
settings.hyper_lambda = 20;
settings.iSTEPminForHyperOptimization = 1;
settings.MaxEvals = '1e6*dim';
settings.MaxEvalsWithSurrogate = '1e4*20';
settings.lambdaMult = 1;
settings.muMult = 1;
settings.largeLambdaMinIter = 3;

settings.withDisp = 0;
settings.maxStepts = 20;
settings.maxerr = 0.45;
settings.alpha = 0.20;
settings.iterstart = 10;
Adapter();
