function run_s_cmaes_bbcomp(FUN, dim, maxfunevals, maxTime)

  persistent archive;

  % FUN = [];
  % dim = [];
  % maxfunevals = [];
  cmParams = struct();
  
  surrogateParams = struct(...
    'evoControl', 'doubletrained', ...
    'evoControlSwitchMode', 'none', ...
    'archive', Archive(dim), ...
    'startTime', tic, ...
    'evoControlMaxTime', maxTime);
  
  xstart = rand(dim, 1);
  
  while maxfunevals > 0
    
    [x, y_evals, stopflag, archive, varargout] = opt_s_cmaes_bbcomp(FUN, dim, maxfunevals, cmParams, surrogateParams, xstart);
    
    maxfunevals = maxfunevals - varargout.evals;
    xstart = cmaesRestartPoint(archive);
    surrogateParams.archive = archive;
  end

end