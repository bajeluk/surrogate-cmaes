function run_s_cmaes_bbcomp(FUN, dim, maxfunevals)

  persistent archive;

  % FUN = [];
  % dim = [];
  % maxfunevals = [];
  cmParams = struct();
  
  sgParams = struct(...
    'evoControl', 'none', ...
    'archive', Archive(dim), ...
    'startTime', tic, ...
    'evoControlMaxTime', 7*24*3600);
  
  nEvals = 0;
  xstart = rand(dim, 1);
  
  while nEvals < maxfunevals
    
    [x, y_evals, stopflag, archive, varargout] = opt_s_cmaes_bbcomp(FUN, dim, maxfunevals, cmParams, sgParams, xstart);
    
    nEvals = size(archive.X, 1);
    
    xstart = cmaesRestartPoint(archive);
  end

end