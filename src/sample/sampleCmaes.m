[fitness_raw] = sampleCmaes(xmean, lambda, BD, N, diagD)

  % Generate and evaluate lambda offspring
 
  fitness.raw = repmat(NaN, 1, lambda + noiseReevals);

  % parallel evaluation
  if flgEvalParallel
      arz = randn(N,lambda);

      if ~flgDiagonalOnly
        arx = repmat(xmean, 1, lambda) + sigma * (BD * arz); % Eq. (1)
      else
        arx = repmat(xmean, 1, lambda) + repmat(sigma * diagD, 1, lambda) .* arz; 
      end

      if noiseHandling 
        if noiseEpsilon == 0
          arx = [arx arx(:,1:noiseReevals)]; 
        elseif flgDiagonalOnly
          arx = [arx arx(:,1:noiseReevals) + ...
                 repmat(noiseEpsilon * sigma * diagD, 1, noiseReevals) ...
                 .* randn(N,noiseReevals)]; 
        else 
          arx = [arx arx(:,1:noiseReevals) + ...
                 noiseEpsilon * sigma * ...
                 (BD * randn(N,noiseReevals))]; 
        end
      end

      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid = arx;
      else
        arxvalid = xintobounds(arx, lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness.raw = feval(fitfun, arxvalid, varargin{:}); 
      countevalNaN = countevalNaN + sum(isnan(fitness.raw));
      counteval = counteval + sum(~isnan(fitness.raw)); 
  end

  % non-parallel evaluation and remaining NaN-values
  % set also the reevaluated solution to NaN
  fitness.raw(lambda + find(isnan(fitness.raw(1:noiseReevals)))) = NaN;  
  for k=find(isnan(fitness.raw)), 
    % fitness.raw(k) = NaN; 
    tries = flgEvalParallel;  % in parallel case this is the first re-trial
    % Resample, until fitness is not NaN
    while isnan(fitness.raw(k))
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
      fitness.raw(k) = feval(fitfun, arxvalid(:,k), varargin{:}); 
      tries = tries + 1;
      if isnan(fitness.raw(k))
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
