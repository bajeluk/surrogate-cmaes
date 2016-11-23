%% Load benchmark 2D data

load 'exp/experiments/exp_doubleEC_11_test/bbob_output/exp_doubleEC_11_test_modellog_8_2D_1.mat'
load 'exp/experiments/exp_doubleEC_11_test/exp_doubleEC_11_test_results_8_2D_1.mat'

g = 10;
last_point_idx = cmaes_out{1}{1}.generationStarts(g+1)-1;
g_generation_idx = [cmaes_out{1}{1}.generationStarts(g):last_point_idx];
orig_idx = logical(cmaes_out{1}{1}.origEvaled(1:last_point_idx));

fgeneric('initialize', 8, 1, '/tmp', struct('algName', 'Test_pure_CMAES'));
model = models{g};

% Preprocess data into GP's familiar matrices

% X_N = cmaes_out{1}{1}.arxvalids(:,find(orig_120));
% y_N = cmaes_out{1}{1}.fvalues(find(orig_120));
% y_N = y_N';

% X_N = (model.trainSigma*model.trainBD) \ (model.dataset.X');
% Trainset in the GP's space
X_N = model.dataset.X';
f_GPToX = @(X) ( model.trainSigma * model.trainBD  * X);
f_XToGP = @(X) ((model.trainSigma * model.trainBD) \ X);

f_yToGP = @(y) (y - model.shiftY) / model.stdY;
f_GPToY = @(y) (y * model.stdY) + model.shiftY;

y_N = f_yToGP(model.dataset.y);
[D, N] = size(X_N);

X_star = f_XToGP(cmaes_out{1}{1}.arxvalids(:,g_generation_idx));
lambda = size(X_star, 2);
[f_star, cov_star] = model.predict(f_GPToX(X_star)');
y_star = fgeneric(f_GPToX(X_star));

mu = floor(lambda/2);

%% Plotting the actual surface (preparation)

landXMin = min(f_GPToX([X_N, X_star]), [], 2);
landXMax = max(f_GPToX([X_N, X_star]), [], 2);
diff = landXMax - landXMin;
landXMin = landXMin - 0.1*diff;
landXMax = landXMax + 0.1*diff;

[landX,landY] = meshgrid(linspace(landXMin(1),landXMax(1),20), linspace(landXMin(2),landXMax(2),20));
landX_l = landX(:); landY_l = landY(:);
landZ_l = zeros(size(landX_l,1),1);
for i = 1:length(landZ_l)
  landZ_l(i) = fgeneric([landX_l(i); landY_l(i)]);
end
landZ = reshape(landZ_l, size(landX,1), size(landX,2));

figure();
surf(landX, landY, landZ);
hold on;
X_N_real = f_GPToX(X_N);
X_star_real = f_GPToX(X_star);
plot3(X_N_real(1,:), X_N_real(2,:), f_GPToY(y_N), 'ro');
plot3(X_star_real(1,:), X_star_real(2,:), y_star, 'ko');

% Plot GP prediction
plot3(X_star_real(1,:), X_star_real(2,:), f_star, 'bo');

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
% % it can be done like this
% Fmu = m_star + K__X_star__X_N * K_inv * (1/sn2) * (y_N - m_N);
% but it is more numerical stable to do it like
Fmu = m_star + K__X_star__X_N * alpha;

% Debug
assert(max(abs(f_GPToY(Fmu) - f_star)) < 1e-6, 'Fmu calculated differs more from model.predict by %e', max(abs(f_GPToY(Fmu) - f_star)));
% fprintf('Fmu calculated differs from model.predict by %e.\n', max(abs(f_GPToY(Fmu) - f_star)));

% [y_GP, std_GP] = gp(model.hyp, model.infFcn, model.meanFcn, model.covFcn, model.likFcn, X_N', y_N, X_star');

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
      M(i,j) = f_GPToY(m_N_p(end) + (h_k - h_l) / (A(j,end) - A(i,end)));
    end
  end
  thresholds = M(:);
  thresholds = sort(thresholds(~isnan(thresholds)));

  %% Identify the rankings in the middle of the inequality boundaries

  % Calculate middle points between the thresholds
  middleThresholds = [2*thresholds(1) - thresholds(end);
    (thresholds(2:end) + thresholds(1:(end-1)))/2;
    2*thresholds(end) - thresholds(1)];

  % Ranking of the most probable Y_s
  mean_rank = ranking(f_GPToY(Fmu(withoutS)))';

  % Determine the ranking before the first threshold
  Y_s = f_yToGP(middleThresholds(1));
  Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
  last_rank = ranking(f_GPToY(Fmu_m))';
  y_ranks = [last_rank];
  rank_diffs = [errRankMu(last_rank, mean_rank, mu)];

  for i = 1:length(thresholds)
    % Find the coordinates of this threshold in the matrix 'M'
    [row, col] = find(M == thresholds(i));
    % This threshold exchanges ranking of this  row <--> col
    this_rank = last_rank;
    this_rank(row) = last_rank(col);
    this_rank(col) = last_rank(row);

    % Save this new ranking and its errRankMu
    y_ranks = [y_ranks; this_rank];
    rank_diffs = [rank_diffs; errRankMu(this_rank, mean_rank, mu)];
    last_rank = this_rank;

    % Debug: just for being sure whether we set the ranking right :)
    Y_s = f_yToGP(middleThresholds(i+1));
    % Calculate vector f* (for this value of Y_s)
    Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; Y_s] - m_N_p);
    % Get the ranking of f*
    check_rank = ranking(f_GPToY(Fmu_m))';
    assert(all(this_rank == check_rank));
  end

  %% Compute the propabilites of these rankings
  % norm_cdfs = [0; normcdf(thresholds, Fmu(s), sqrt(Fs2(s)/max(Fs2))); 1];
  norm_cdfs = [0; normcdf(thresholds, f_GPToY(Fmu(s)), Fs2(s)*model.stdY); 1];
  probs = norm_cdfs(2:end) - norm_cdfs(1:(end-1));

  % Merge expected errors for corresponding rankings
  rankings_errs = containers.Map();
  for r = 1:length(probs)
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

  %% Plot the two extreme rankings for chosen s
  X_star_m_real = f_GPToX(X_star_m);
  X_s_real = f_GPToX(X_s);

  if (s == 6)
    y_extremes = [thresholds(1)-sqrt(eps), thresholds(end)+sqrt(eps)];
    colors = ['m', 'g'];
    for yi = 1:2
      y = f_yToGP(y_extremes(yi));
      % Ks = feval(model.covFcn{:}, model.hyp.cov, X_N_p', x_surf')';
      % ms_extreme = feval(model.meanFcn, model.hyp.mean, x_surf');
      % y_extreme = ms_extreme + Ks * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p);
      % plot(x_surf, y_extreme, [colors(yi) '--']);
      Fmu_m = m_star_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N; y] - m_N_p);
      plot3(X_star_m_real(1,:), X_star_m_real(2,:), f_GPToY(Fmu_m), [colors(yi) '+'])
      plot3(X_s_real(1,:), X_s_real(2,:), y_extremes(yi), 'r*');
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
