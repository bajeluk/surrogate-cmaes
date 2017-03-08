function [perm, errs] = expectedRankDiff(model, arxvalid, mu, varargin)
%EXPECTEDRANKDIFF returns permutation of points which causes highest expected difference in ranking
%
% model   the best Gaussian process model from current generation
% arxvalid  population of new points from which should be chosen
% mu      CMA-ES' mu parameter -- how many points are taken for
%         updating mean and rank-mu update of covariance matrix
%
% perm    Returns permutation of the points in the order from maximum
%         expected error (in perm(1)) to the lowest error (in
%         perm(end)
% errs    Expected errors themselves (in the order corresponding to arxvalid)

% TODO:
% [ ] rewrite the core calculation of ranking differences into the MEX file

  % if (nargin >= 5)
  %   rankFunc = varargin{1};
  % else
  %   rankFunc = @errRankMu;
  % end
  rankFunc = @errRankMu;

  % conversion functions from/to space of GP model and real space
  f_GPToX = @(X) ( model.trainSigma * model.trainBD  * X);
  f_XToGP = @(X) ((model.trainSigma * model.trainBD) \ X);
  f_yToGP = @(y) (y - model.shiftY) / model.stdY;
  f_GPToY = @(y) (y * model.stdY) + model.shiftY;

  % matrix of 'training' data points (already converted to model-space)
  X_N = model.getDataset_X()';
  % vector of 'training' f-values
  % (model.getDataset_y() is not yet converted to model-space)
  y_N = f_yToGP(model.getDataset_y());
  % number of 'training' datapoints
  N = size(X_N, 2);

  % population of new points in the model-space
  X_star = f_XToGP(arxvalid);
  % population size
  lambda = size(X_star, 2);

  % Debug: for cross-check purposes whether we have done the GP prediction right
  [f_star, cov_star] = model.predict(arxvalid');

  % Calculate covariances
  %
  K__X_star__X_N = feval(model.covFcn{:}, model.hyp.cov, X_N', X_star')'; % cross-covariances
  K__X_N__X_N = feval(model.covFcn{:}, model.hyp.cov, X_N');    % train covariance (posterior?)
  Kss = feval(model.covFcn{:}, model.hyp.cov, X_star', 'diag'); % self-variance

  % Posterior
  %
  % evaluate mean vector for X_N
  m_N = feval(model.meanFcn, model.hyp.mean, X_N');
  % noise variance of likGauss
  sn2 = exp(2*model.hyp.lik);
  % Cholesky factor of covariance with noise
  L = chol(K__X_N__X_N/sn2 + eye(N) + 0.0001*eye(N));
  % inv(K+noise) * (y_N - mean)
  alpha = solve_chol(L, y_N - m_N) / sn2;

  % Predictive mean
  %
  % evaluate mean function for test points X_star
  m_star = feval(model.meanFcn, model.hyp.mean, X_star');
  % Predictive conditional mean fs|f (shorter form, using calculated alpha)
  Fmu = m_star + K__X_star__X_N * alpha;
  % Predictive conditional mean fs|f (longer form which uses K_inv)
  % K_inv = solve_chol(L, eye(N));        % Inverse of the covariance
  % Fmu = m_star + K__X_star__X_N * K_inv * (1/sn2) * (y_N - m_N);

  % Debug
  % assert(abs(max(f_GPToY(Fmu) - f_star)/max(f_star)) < 1e-6, 'Fmu calculated relatively differs from model.predict by factor %e', abs(max(f_GPToY(Fmu) - f_star)/max(f_star)));
  if (abs(max(f_GPToY(Fmu) - f_star)/max(f_star)) > 1e-6)
    fprintf(2, 'Fmu calculated relatively differs from model.predict by factor %e', abs(max(f_GPToY(Fmu) - f_star)/max(f_star)));
  end

  % Predictive variances
  %
  Ltril = all(all(tril(L,-1)==0));   % is L an upper triangular matrix?
  if (Ltril)   % L is triangular => use Cholesky parameters (alpha,sW,L)
    sW = ones(N,1) / sqrt(sn2);      % sqrt of noise precision vector
    V  = L'\(repmat(sW,1,lambda) .* K__X_star__X_N');
    fs2 = Kss - sum(V.*V,1)';        % predictive variances
  else         % L is not triangular => use alternative parametrisation
    fs2 = Kss + sum(K__X_star__X_N .* (L*K__X_star__X_N'), 1)';  % predictive variances
  end
  % remove numerical noise i.e. negative variances
  Fs2 = max(fs2, 0);
  % apply likelihood function
  [~, ~, Ys2] = feval(model.likFcn, model.hyp.lik, [], Fmu, Fs2);
  % correct GPML's bug for the case when Fs2 == zeros(...)
  if (size(Ys2,1) == 1)
    Ys2 = repmat(Ys2, size(Fs2,1), 1);
  end

  % Debug
  % assert(abs(max(Ys2*model.stdY^2 - cov_star)/max(cov_star)) < 1e-4, 'Ys2 calculated differs more from model.predict by %e \%', abs(max(Ys2*model.stdY^2 - cov_star)/max(cov_star)));
  if (abs(max(Ys2*(model.stdY^2) - cov_star)/max(cov_star)) > 1e-6)
    fprintf(2, 'Ys2 calculated relatively differs from model.predict by factor %e\n', abs(max(Ys2*model.stdY^2 - cov_star)/max(cov_star)));
  end

  % Iterate through all lambda points
  %
  expectedErr = zeros(lambda,1);

  for s = 1:lambda

    % Calculate modified covariances and prediction as if 'x_s' would be in the train set
    withoutS = [1:(s-1), (s+1):lambda]';
    X_star_m = X_star(:,withoutS);
    X_s      = X_star(:,s);
    X_N_p    = [X_N X_s];

    % covariances
    K__X_star_m__X_N_p = feval(model.covFcn{:}, model.hyp.cov, X_N_p', X_star_m')';
    new_vect = feval(model.covFcn{:}, model.hyp.cov, X_s', X_N');
    K__X_N_p__X_N_p = [K__X_N__X_N, new_vect'; new_vect, Kss(end)];
    % the two previous lines are the same as this:
    % K__X_N_p__X_N_p = feval(model.covFcn{:}, model.hyp.cov, X_N_p', []);
    % evaluate mean vectors
    m_N_p  = [m_N; m_star(s)];
    % m_N_p  = feval(model.meanFcn, model.hyp.mean, X_N_p');
    m_star_m = m_star(withoutS);
    % m_star_m = feval(model.meanFcn, model.hyp.mean, X_star_m');
    % noise variance of likGauss
    % already calculated: sn2 = exp(2*model.hyp.lik);
    % Cholesky factor of covariance with noise
    % [Lp, p] = chol(K__X_N_p__X_N_p/sn2+eye(N+1) + 0.0001*eye(N+1));
    [Lp, p] = chol(K__X_N_p__X_N_p/sn2);
    if (p > 0)
      expectedErr(s) = 0;
      fprintf(2, 'expectedRankDiff(): GP''s K__X_N_p covariance matrix is not positive definite for inverse. Setting err(%d) = 0.\n', s);
      continue;
    end

    % Covariance * y_N
    Kp_inv = solve_chol(Lp, eye(N+1));

    % Calculate inequalitiy boudnaries
    %
    A = K__X_star_m__X_N_p * Kp_inv * (1/sn2);
    M = NaN(lambda-1,lambda-1);
    M2 = NaN(lambda-1,lambda-1);
    h = - m_star_m + A(:,1:(end-1)) * (y_N - m_N);
    for i = 1:(lambda-1)
      for j = (i+1):(lambda-1)
        M(i,j) = m_N_p(end) + (h(i) - h(j)) / (A(j,end) - A(i,end));
      end
    end
    M = f_GPToY(M);
    thresholds = M(:);
    [thresholds, tidx] = sort(thresholds);

    nThresholds = (lambda-1)*(lambda-2)/2;
    thresholds = thresholds(1:nThresholds);
    tidx       = tidx(1:nThresholds);

    % Calculate middle points between the thresholds
    middleThresholds = [2*thresholds(1) - thresholds(end);
      (thresholds(2:end) + thresholds(1:(end-1)))/2;
      2*thresholds(end) - thresholds(1)];

    % Determine the different rankings in each interval between thresholds
    %
    % Ranking of the most probable Y_s == f*(s)
    mean_rank = ranking(f_GPToY(Fmu(withoutS)))';

    % Determine the first ranking, i.e. before the first threshold
    Y_s = f_yToGP(middleThresholds(1));
    Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
    this_rank = ranking(f_GPToY(Fmu_m))';
    y_ranks      = zeros(nThresholds, lambda-1);
    y_ranks(1,:) = this_rank;
    rank_diffs    = zeros(nThresholds+1, 1);
    mu = min(mu, lambda-1);
    [rank_diffs(1), ~, maxErr] = rankFunc(this_rank, mean_rank, mu);
    % This is for reducing complexity of errRankMu error calculations
    [~, si_mean_rank] = sort(mean_rank);
    % rows = mod(tidx, lambda-1);               % not significant speed-up
    % cols = floor(tidx/(lambda-1)) + 1;        % not significant speed-up

    % TODO: rewrite this as a MEX function
    %
    % function calculateRankDiffs(lambda, mu, tidx, this_rank, mean_rank)
    %
    % This is the critical section -- the inner of the for-cycle
    % is called O(lambda^3) per generation!
    % -- for each of the point (cycle through 's') there are O(lambda^2) thresholds
    %
    for i = 1:nThresholds
      % Find the coordinates of the i-th threshold in the matrix 'M'
      row = mod(tidx(i), lambda-1);
      col = floor(tidx(i)/(lambda-1)) + 1;
      % This threshold exchanges ranking of this  row <--> col, do it in ranking
      tmp = this_rank(row);
      this_rank(row) = this_rank(col);
      this_rank(col) = tmp;

      % Save this new ranking and its errRankMu
      y_ranks(i+1,:)  = this_rank;

      % Calculate the new errRankMu error
      % This is computationally expensive -- O(lambda * log(lambda))
      %   rank_diffs(i+1) = rankFunc(this_rank, mean_rank, mu);

      % This is cheaper -- O(mu) -- but calling a function still costs a lot!
      %   rank_diffs(i+1) = errRankMuLite(si_this_rank, si_mean_rank, mu, maxErr);

      % And this is even cheaper, without calling a function -- O(mu)
      % -- for lambda=35 two times faster!
      %
      % take the first 'mu' elements of 'this_ranking' sorted according to
      % the mean_rank's ranking
      r1 = this_rank(si_mean_rank(1:mu));
      % calculate the difference to the ranking 1:mu
      % and leave the normalization after the for-cycle
      rank_diffs(i+1) = sum(abs((1:mu) - r1));

      % Debug:
      % assert(rank_diffs(i+1)/maxErr == rankFunc(this_rank, mean_rank, mu), 'expectedRankDiff -- errRank is wrong!');
      % assert(rank_diffs(i+1) == rankFunc(this_rank, mean_rank, mu), 'expectedRankDiff -- errRank is wrong!');

      % Debug: just for being sure whether we set the ranking right :)
      % Y_s = f_yToGP(middleThresholds(i+1));
      % % Calculate vector f* (for this value of Y_s)
      % Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
      % % Get the ranking of f*
      % check_rank = ranking(f_GPToY(Fmu_m))';
      % assert(all(this_rank == check_rank), 'expectedRankDiff -- ranking is wrong!');
    end
    rank_diffs(2:end) = rank_diffs(2:end) ./ maxErr;

    % Compute the propabilites of these rankings
    %
    % norm_cdfs = [0; normcdf(thresholds, Fmu(s), sqrt(Ys2(s)/max(Ys2))); 1];
    norm_cdfs = [0; normcdf(thresholds, f_GPToY(Fmu(s)), Ys2(s)*model.stdY^2); 1];
    probs = norm_cdfs(2:end) - norm_cdfs(1:(end-1));

    % Save the resulting expected error for this individual 's'
    expectedErr(s) = sum(rank_diffs.*probs);

  end  % for s = 1:lambda

  % Return the permutation with the highest expected 
  % rankDiff error of points according to the descending errors
  [~, perm] = sort(-expectedErr);
  errs = expectedErr;

  [~, sd2i] = sort(-Ys2);
  nPoints = 2;
  if (all(sd2i(1:nPoints) == perm(1:nPoints)))
    fprintf(2, 'The first %d points in "expectedRankDiff" are the same as in the "sd2" criterion.\n', nPoints);
  end
end
