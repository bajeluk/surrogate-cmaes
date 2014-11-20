function [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaesOnlyFitness(arx, arxvalid, arz, xmean, sigma, lambda, BD, diagD, fitfun, fitargs, opts)
% sampleCmaesOnlyFitness  evaluation of the individuals in @arx/@arxvalid with fitness
%
% It generates new samples for the individuals for which NaN function value is
% returned according to @xmean, @sigma, @BD and @diagD

  % CMA-ES sampling options
  noiseReevals = opts.noiseReevals;
  bnd.isactive = opts.isBoundActive;
  lbounds = opts.lbounds;
  ubounds = opts.ubounds;
  counteval = opts.counteval;
  flgEvalParallel = opts.flgEvalParallel;
  flgDiagonalOnly = opts.flgDiagonalOnly;
  xintobounds = opts.xintobounds;

  isBoundActive = any(lbounds > -Inf) || any(ubounds < Inf); 
  N = size(xmean, 1);

  % Generate and evaluate lambda offspring
 
  fitness_raw = NaN(1, lambda + noiseReevals);
  countevalNaN = 0;

  % parallel evaluation
  if flgEvalParallel

      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness_raw = feval(fitfun, arxvalid, fitargs{:}); 
      countevalNaN = countevalNaN + sum(isnan(fitness_raw));
      counteval = counteval + sum(~isnan(fitness_raw)); 
  end

  % non-parallel evaluation and remaining NaN-values
  % set also the reevaluated solution to NaN
  fitness_raw(lambda + find(isnan(fitness_raw(1:noiseReevals)))) = NaN;  
  for k=find(isnan(fitness_raw)), 
    % fitness_raw(k) = NaN; 
    tries = flgEvalParallel;  % in parallel case this is the first re-trial
    % Resample, until fitness is not NaN
    while isnan(fitness_raw(k))
      if k <= lambda  % regular samples (not the re-evaluation-samples)
        arz(:,k) = randn(N,1); % (re)sample

        if flgDiagonalOnly  
          arx(:,k) = xmean + sigma * diagD .* arz(:,k);              % Eq. (1)
        else
          arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)
        end
      else % re-evaluation solution with index > lambda
        if flgDiagonalOnly  
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * diagD .* randn(N,1);
        else
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * (BD * randn(N,1));
        end
      end
      
      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid(:,k) = arx(:,k);
      else
        arxvalid(:,k) = xintobounds(arx(:,k), lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx should not be changed.
      fitness_raw(k) = feval(fitfun, arxvalid(:,k), fitargs{:}); 
      tries = tries + 1;
      if isnan(fitness_raw(k))
	countevalNaN = countevalNaN + 1;
      end
      if mod(tries, 100) == 0
	warning([num2str(tries) ...
                 ' NaN objective function values at evaluation ' ...
                 num2str(counteval)]);
      end
    end
    counteval = counteval + 1; % retries due to NaN are not counted
  end
end
