function exp_saACMES_task(taskpath, funs_, varargin)

  cd([taskpath filesep '..' filesep '..']);
  startup;
  cd([taskpath filesep 'saACMESlambdaRevMinIter3v2']);
  disp('Current path:');
  pwd

  global settings;

  if (~isempty(varargin))
    dims_ = varargin{1};
  else
    dims_ = [2, 3, 5, 10, 20];
  end

  fprintf(1, 'Functions : %d\n', funs_);
  fprintf(1, 'Dimensions: %d\n', dims_);
settings.instances = [1:5,31:40];
  settings.dims = dims_;
  settings.funs = funs_;
  settings.pathname = 'd1';
  settings.algname = '';
  settings.ntarray = [1];
  settings.savfile = 'r1';

  settings.BIPOP = 0; 
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
  settings.MaxEvals = '250*dim';
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
end
