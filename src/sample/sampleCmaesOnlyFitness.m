function [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaesOnlyFitness(arx, arxvalid, arz, sigma, lambda, counteval, cmaesState, opts, varargin)
% sampleCmaesOnlyFitness  evaluation of the individuals in @arx/@arxvalid with fitness
%
% It generates new samples for the individuals for which NaN function value is
% returned according to @xmean, @sigma, @BD and @diagD

  MAX_TRIES = 1000;

  % CMA-ES state variables
  xmean = cmaesState.xmean;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  fitfun= cmaesState.fitfun_handle;

  % CMA-ES sampling options
  noiseReevals = opts.noiseReevals;
  bnd.isactive = opts.isBoundActive;
  lbounds = opts.lbounds;
  ubounds = opts.ubounds;
  flgEvalParallel = opts.flgEvalParallel;
  flgDiagonalOnly = opts.flgDiagonalOnly;
  xintobounds = opts.xintobounds;

  % input options
  archive = [];
  if (nargin >= 6 && ~isempty(varargin) && ischar(varargin{1}))
    switch lower(varargin{1})
      case 'archive'
        archive = varargin{2};
        varargin(1:2) = [];
    end
  end

  % cut lambda to the maximal number of allowed function evaluations
  orig_lambda = lambda;
  maxevals = defopts(cmaesState, 'thisGenerationMaxevals', lambda + noiseReevals);
  maxevals = min(maxevals, lambda+noiseReevals);
  lambda = min(maxevals, lambda);
  noiseReevals = min(noiseReevals, maxevals - lambda);

  isBoundActive = any(lbounds > -Inf) || any(ubounds < Inf);
  N = size(xmean, 1);

  % Generate and evaluate lambda offspring

  fitness_raw = NaN(1, orig_lambda + noiseReevals);
  countevalNaN = 0;

  isEvaled = false(1, orig_lambda + noiseReevals);
  if (~isempty(archive))
    for k = 1:size(fitness_raw, 2)
      % Do not re-evaluate the already saved point, but use
      % the value from archive
      [isAlreadySaved, idx] = archive.isInArchive(arxvalid(:,k)');
      if (isAlreadySaved)
        fitness_raw(k) = archive.y(idx);
        isEvaled(k) = true;
      end
    end
  end

  % do not evaluate more than lambda+noiseReevals points
  toOrigEval = ~isEvaled;
  toOrigEval((lambda+noiseReevals+1):end) = false;

  % parallel evaluation
  if flgEvalParallel

      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness_raw(toOrigEval) = feval(fitfun, arxvalid(toOrigEval), varargin{:});
  else
      for k = 1:size(fitness_raw, 2)
        if (toOrigEval(k))
          fitness_raw(k) = feval(fitfun, arxvalid(:,k), varargin{:});
        end
      end
  end
  countevalNaN = countevalNaN + sum(isnan(fitness_raw(toOrigEval)));
  counteval = counteval + sum(~isnan(fitness_raw(toOrigEval)));

  % non-parallel evaluation and remaining NaN-values
  % set also the reevaluated solution to NaN
  fitness_raw(lambda + find(isnan(fitness_raw(1:noiseReevals)))) = NaN;
  nanPoints = find(isnan(fitness_raw(1:(orig_lambda+noiseReevals))));
  for k = nanPoints
    % fitness_raw(k) = NaN;
    tries = 1;  % we have already tried f-eval once (see above)
    % Resample, until fitness is not NaN
    while (isnan(fitness_raw(k)) && tries < MAX_TRIES)
      if k <= orig_lambda  % regular samples (not the re-evaluation-samples)
        arz(:,k) = randn(N,1); % (re)sample

        if flgDiagonalOnly
          arx(:,k) = xmean + sigma * diagD .* arz(:,k);              % Eq. (1)
        else
          arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)
        end
      else % re-evaluation solution with index > lambda
        if flgDiagonalOnly
          arx(:,k) = arx(:,k-orig_lambda) + (noiseEpsilon * sigma) * diagD .* randn(N,1);
        else
          arx(:,k) = arx(:,k-orig_lambda) + (noiseEpsilon * sigma) * (BD * randn(N,1));
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
      if (~isempty(archive))
        [isAlreadySaved, idx] = archive.isInArchive(arxvalid(:,k)');
        if (isAlreadySaved)
          fitness_raw(k) = archive.y(idx);
          counteval = counteval - 1;
          break;
        end
      end

      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx should not be changed.

      % orignal-evaluate only if you did not exceed 'lambda'
      if (k <= lambda + noiseReevals)
        fitness_raw(k) = feval(fitfun, arxvalid(:,k), varargin{:});
        tries = tries + 1;
        if isnan(fitness_raw(k))
          countevalNaN = countevalNaN + 1;
        end
        if mod(tries, 100) == 0
          warning([num2str(tries) ...
                  ' NaN objective function values at evaluation ' ...
                  num2str(counteval)]);
        end
      else
        % Mock fitness for points after the limit of MaxEvals
        fitness_raw(k) = Inf;
        counteval = counteval - 1;
        break;
      end
    end
    counteval = counteval + 1; % retries due to NaN are not counted
  end
end
