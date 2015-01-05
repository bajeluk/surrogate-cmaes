function xbest = MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% samples new points uniformly randomly in [-5,5]^DIM
% and evaluates them on FUN until ftarget of maxfunevals
% is reached, or until 1e8 * DIM fevals are conducted. 

  maxfunevals = min(1e8 * DIM, maxfunevals); 
  popsize = min(maxfunevals, 200);
  fbest = inf;
  for iter = 1:ceil(maxfunevals/popsize)
    xpop = 10 * rand(DIM, popsize) - 5;      % new solutions
    [fvalues, idx] = sort(feval(FUN, xpop)); % evaluate
    if fbest > fvalues(1)                    % keep best
      fbest = fvalues(1);
      xbest = xpop(:,1);
    end
    if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
      break;                                 % (works also for noisy functions)
    end
  end 

  
