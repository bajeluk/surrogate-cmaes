%% Generate a benchmark 1D data

if (exist('rankdiff_crit_test_1D.mat', 'file'))
  load('rankdiff_crit_test_1D.mat', 'f', 'trainSize', 'lambda', 'mu', 'X_N', 'y_N', 'D', 'N', 'X_star');
else
  f = @(x) 0.25*(x+2).*sin(pi*x).^4;

  trainSize = 7;
  lambda = 6;
  mu = 3;

  X_N = -0.2 + randn(1,trainSize)*0.4;
  y_N = f(X_N)';
  % set some noise on the training data
  y_N = y_N +  randn(size(y_N))*0.02;
  [D, N] = size(X_N);

  X_star = 0.2 + randn(1,lambda)*0.2;
  lambda = size(X_star, 2);
end

%% Plot initial settings

% sort the population from left to right
X_star = sort(X_star);

x_surf = -1:0.02:1;
figure();
plot(x_surf, f(x_surf), 'k-');
ax = gca;
ax.YLim = [-0.5 1];
hold on;
plot(X_N, y_N, 'r+');
plot(X_star, -0.5*ones(size(X_star)), 'b*');

%% Construct GP model

modelOpts = struct('transformCoordinates', false, 'normalizeY', false);
model = GpModel(modelOpts, 0);
cmState = struct('xmean', -0.2, 'countiter', 1, 'sigma', 0.3, 'lambda', lambda, 'BD', 1, 'diagD', 1);
sampleOpts = struct('noiseReevals', 0, 'isBoundActive', 1, 'lbounds', -1, 'ubounds', 1, ...
  'flgDiagonalOnly', false, 'noiseHandling', false, 'xintobounds', @xintobounds);
model = model.train(X_N', y_N, cmState, sampleOpts);
model.trainGeneration = 1;
if (~iscell(model.likFcn))
  model.likFcn = {model.likFcn};
end
y_star = f(X_star');
[f_star, cov_star] = model.predict(X_star');

%% Plot GP prediction

plot(X_star, f_star, 'b+');
y_model = model.predict(x_surf');
plot(x_surf, y_model', 'b--');

%% Evaluate covariance function

K__X_star__X_N = feval(model.covFcn{:}, model.hyp.cov, X_N', X_star')'; % cross-covariances
K__X_N__X_N = feval(model.covFcn{:}, model.hyp.cov, X_N');    % train covariance (posterior?)

Kss = feval(model.covFcn{:}, model.hyp.cov, X_star', 'diag'); % self-variance

%% Posterior

% evaluate mean vector for X_N
m_N = feval(model.meanFcn, model.hyp.mean, X_N');
% noise variance of likGauss
sn2 = exp(2*model.hyp.lik);
% Cholesky factor of covariance with noise
L = chol(K__X_N__X_N/sn2 + eye(N) + 0.0001*eye(N));
% inv(K+noise) * (y_N - mean)
alpha = solve_chol(L, y_N - m_N) / sn2;

%% Predictive mean

% evaluate mean function for test points X_star
m_star = feval(model.meanFcn, model.hyp.mean, X_star');

% conditional mean fs|f (shorter form)
% Fmu = ms + K__X_star__X_N * alpha;

% Inverse of the covariance
K_inv = solve_chol(L, eye(N));

% Predictive conditional mean fs|f
Fmu = m_star + K__X_star__X_N * K_inv * (1/sn2) * (y_N - m_N);

%% Predictive variances

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

%% Iterate through all lambda points
expectedErr = zeros(lambda,1);

for s = 1:lambda
  %% Calculate modified covariances and prediction as if 'x_s' would be in the train set
  withoutS = [1:(s-1), (s+1):lambda]';

  X_star_m = X_star(:,withoutS);
  X_s      = X_star(:,s);
  X_N_p    = [X_N X_s];
  % plot(X_s, Fmu(s), 'b*');

  % covariances
  K__X_star_m__X_N_p = feval(model.covFcn{:}, model.hyp.cov, X_N_p', X_star_m')';
  K__X_N_p__X_N_p = feval(model.covFcn{:}, model.hyp.cov, X_N_p', []);
  % evaluate mean vectors
  m_N_p  = feval(model.meanFcn, model.hyp.mean, X_N_p');
  m_star_m = feval(model.meanFcn, model.hyp.mean, X_star_m');
  % noise variance of likGauss
  sn2 = exp(2*model.hyp.lik);
  % Cholesky factor of covariance with noise
  Lp = chol(K__X_N_p__X_N_p/sn2+eye(N+1) + 0.0001*eye(N+1));
  % Covariance * y_N
  Kp_inv = solve_chol(Lp, eye(N+1));

  %% Calculate inequalitiy boudnaries

  A = K__X_star_m__X_N_p * Kp_inv * (1/sn2);

  M = NaN(lambda-1,lambda-1);

  for i = 1:(lambda-1)
    for j = (i+1):(lambda-1)
      h_k = - m_star_m(i) + A(i,1:(end-1)) * (y_N - m_N);
      h_l = - m_star_m(j) + A(j,1:(end-1)) * (y_N - m_N);
      M(i,j) = m_N_p(end) + (h_k - h_l) / (A(j,end) - A(i,end));
    end
  end
  thresholds = M(:);
  thresholds = sort(thresholds(~isnan(thresholds)));

  %% Identify the rankings in the middle of the inequality boundaries
  
  middleThresholds = [2*thresholds(1) - thresholds(end);
    thresholds(2:end) - thresholds(1:(end-1));
    2*thresholds(end) - thresholds(1)];
  
  y_ths   = [];
  y_ranks = [];
  rank_diffs = [];
  mean_rank = ranking(Fmu(withoutS))';
  for i = 1:length(thresholds)
    Y_s = middleThresholds(i);
    % Calculate vector f* (for this value of Y_s)
    Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
    % Get the ranking of f*
    this_rank = ranking(Fmu_m)';
    y_ranks = [y_ranks; this_rank];
    y_ths = [y_ths; thresholds(i)];
    rank_diffs = [rank_diffs; errRankMu(this_rank, mean_rank, mu)];
  end
  % Save also the last ranking
  Y_s = middleThresholds(end);
  Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
  y_ranks = [y_ranks; ranking(Fmu_m)'];
  rank_diffs = [rank_diffs; errRankMu(this_rank, mean_rank, mu)];

  % %% Sample possible values of Y_s and identify different orderings of f*
  % 
  % % Try different values of Y_s which has prior
  % %   Y_s ~ N(Fmu(s), Fs2(s))
  % 
  % n_sample = 5000;
  % sigma_factor = 20;
  % 
  % lb = Fmu(s) - sigma_factor*Fs2(s);
  % ub = Fmu(s) + sigma_factor*Fs2(s);
  % y_space = linspace(lb, ub, n_sample);
  % mean_rank = ranking(Fmu(withoutS))';
  % y = lb;
  % this_rank = ranking(m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p))';
  % last_rank = this_rank;
  % last_y  = lb;
  % y_ths   = [];
  % y_ranks = [];
  % rank_diffs = [];
  % 
  % % Iterate through different values of Y_s
  % %   y_ths is the last value of the interval where y_ranks holds
  % for y = y_space(2:end)
  %   % Calculate vector f* (for this value of Y_s)
  %   Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p);
  %   % Get the ranking of f*
  %   this_rank = ranking(Fmu_m)';
  %   if (~all(this_rank == last_rank))
  %     % The ranking has changed from the last iteration => save it
  %     y_ths = [y_ths; last_y];
  %     y_ranks = [y_ranks; last_rank];
  %     rank_diffs = [rank_diffs; errRankMu(last_rank, mean_rank, mu)];
  %     last_rank = this_rank;
  %   end
  %   last_y    = y;
  % end
  % % Save also the last ranking
  % y_ranks = [y_ranks; this_rank];
  % rank_diffs = [rank_diffs; errRankMu(this_rank, mean_rank, mu)];

  %% Compute the propabilites of these rankings
  norm_cdfs = [0; normcdf(y_ths, Fmu(s), Fs2(s)); 1];
  probs = norm_cdfs(2:end) - norm_cdfs(1:(end-1));
  
  % Merge expected errors for corresponding rankings
  rankings_errs = containers.Map();
  for r = 1:length(probs)
    if (probs(r) < 1e-10)
      continue;
    end
    key = sprintf('%d ', y_ranks(r,:));
    if (rankings_errs.isKey(key))
      rankings_errs(key) = rankings_errs(key) + rank_diffs(r)*probs(r);
    else
      rankings_errs(key) = rank_diffs(r)*probs(r);
    end
  end
  errors = rankings_errs.values;
  delete(rankings_errs);
  
  % Save the resulting expected error for this individual 's'
  expectedErr(s) = sum([errors{:}]);
  
  %% Plot the two extreme rankings for chosen s == 4
  if (s == 4)
    y_extremes = [y_ths(1), y_ths(end)+0.02];
    colors = ['m', 'g'];
    for yi = 1:2
      y = y_extremes(yi);
      Ks = feval(model.covFcn{:}, model.hyp.cov, X_N_p', x_surf')';
      ms_extreme = feval(model.meanFcn, model.hyp.mean, x_surf');
      y_extreme = ms_extreme + Ks * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p);
      plot(x_surf, y_extreme, [colors(yi) '--']);
      Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p);
      plot(X_star_m, Fmu_m, [colors(yi) '+'])
      plot(X_s, y, 'r*');
    end

    legend('True function', 'GP train set of size N (with small noise)', '\lambda Sampled points (on x-axis)', ...
      '\lambda Sampled points (mean predicted)', 'Mean prediction based on the N-train set', ...
      'Extreme ordering GP with (N+1)-trainset');
  end
end

disp('Expected errors for the possible choices of Y_s:');
disp(expectedErr);
 
% %% Topological sort
% 
% nt = size(thresholds,1);
% S = -1 * ones(nt, lambda-1);
% Ok = false(nt,1);
% for t = 1:nt
%   G = (M < thresholds(t)) + (M' >= thresholds(t));
%   [S(t,:), Ok(t)] = toposort(G);
%   
%   Y_s = thresholds(t)+sqrt(eps);
%   Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
%   disp([Ok(t) ranking(Fmu_m)'; NaN S(t,:)]);
% end
