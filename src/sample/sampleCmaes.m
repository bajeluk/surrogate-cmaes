function [fitness_raw, arx, arxvalid, arz, counteval_out] = sampleCmaes(cmaesState, opts, lambda, counteval_in, varargin)

  % TODO: rewrite the meaning of counteval as the number of _NEW_ original
  %       evaluations made during this specific call

  % input options
  if (nargin >= 6 && ischar(varargin{1}) && strcmpi(varargin{1}, 'Archive'))
    archive = varargin{2};
    varargin(1:2) = [];
  else
    archive = [];
  end

  % CMA-ES state variables
  xmean = cmaesState.xmean;
  sigma = cmaesState.sigma;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  fitfun = cmaesState.fitfun_handle;
  % Note: lambda is a variable parameter defining # of point to sample,
  % it is not always:
  %   lambda == cmaesState.lambda

  % CMA-ES sampling options
  noiseReevals = opts.noiseReevals;
  bnd.isactive = opts.isBoundActive;
  lbounds = opts.lbounds;
  ubounds = opts.ubounds;
  flgEvalParallel = opts.flgEvalParallel;
  flgDiagonalOnly = opts.flgDiagonalOnly;
  noiseHandling = opts.noiseHandling;
  xintobounds = opts.xintobounds;

  isBoundActive = any(lbounds > -Inf) || any(ubounds < Inf); 
  N = size(xmean, 1);

  % Generate and evaluate lambda offspring

  fitness_raw = NaN(1, lambda + noiseReevals);
  countevalNaN = 0;
  newCounteval = 0;

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

      if ~isBoundActive
        arxvalid = arx;
      else
        arxvalid = xintobounds(arx, lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness_raw = feval(fitfun, arxvalid, varargin{:});
      countevalNaN = countevalNaN + sum(isnan(fitness_raw));
      newCounteval = sum(~isnan(fitness_raw));
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

      % Do not re-evaluate the already saved point, but use
      % the value from archive
      if (~isempty(archive) && k <= lambda)
        [isAlreadySaved, idx] = archive.isInArchive(arxvalid(:,k)');
        if (isAlreadySaved)
          fitness_raw(k) = archive.y(idx);
          newCounteval = newCounteval - 1;
          break;
        end
      end

      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx should not be changed.
      fitness_raw(k) = feval(fitfun, arxvalid(:,k), varargin{:});
      tries = tries + 1;
      if isnan(fitness_raw(k))
	      countevalNaN = countevalNaN + 1;
      end
      if mod(tries, 100) == 0
	      warning([num2str(tries) ...
                 ' NaN objective function values at evaluation ' ...
                 num2str(counteval_in + newCounteval)]);
      end
    end
    newCounteval = newCounteval + 1; % retries due to NaN are not counted
  end

  counteval_out = counteval_in + newCounteval;
end
