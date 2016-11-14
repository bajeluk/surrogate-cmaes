%% Load data saved from DTS-CMA-ES run

load 'exp/experiments/exp_doubleEC_11_test/bbob_output/exp_doubleEC_11_test_modellog_8_2D_1.mat'
load 'exp/experiments/exp_doubleEC_11_test/exp_doubleEC_11_test_results_8_2D_1.mat'

g = 10;
last_point_idx = cmaes_out{1}{1}.generationStarts(g+1)-1;
g_generation_idx = [cmaes_out{1}{1}.generationStarts(g):last_point_idx];
orig_idx = logical(cmaes_out{1}{1}.origEvaled(1:last_point_idx));

fgeneric('initialize', 8, 1, '/tmp', struct('algName', 'Test_pure_CMAES'));

%% Preprocess data into GP's familiar matrices

% X_N = cmaes_out{1}{1}.arxvalids(:,find(orig_120));
% y_N = cmaes_out{1}{1}.fvalues(find(orig_120));
% y_N = y_N';

X_N = models{g}.trainSigma*models{g}.trainBD*models{g}.dataset.X';
y_N = models{g}.dataset.y;
[D N] = size(X_N);

X_star = cmaes_out{1}{1}.arxvalids(:,g_generation_idx);
lambda = size(X_star, 2);
[f_star, cov_star] = models{g}.predict(X_star');
y_star = fgeneric(X_star);

%% Plotting the actual surface (preparation)

landXMin = min([X_N, X_star], [], 2);
landXMax = max([X_N, X_star], [], 2);
diff = landXMax - landXMin;
landXMin = landXMin - 0.1*diff;
landXMax = landXMax + 0.1*diff;

[landX,landY] = meshgrid(landXMin(1):0.2:landXMax(1), landXMin(2):0.2:landXMax(2));
landX_l = landX(:); landY_l = landY(:);
landZ_l = zeros(size(landX_l,1),1);
for i = 1:length(landZ_l)
  landZ_l(i) = fgeneric([landX_l(i); landY_l(i)]);
end
landZ = reshape(landZ_l, size(landX,1), size(landX,2));

%% Evaluate covariance function

K__X_star__X_N = feval(models{g}.covFcn{:}, models{g}.hyp.cov, X_N', X_star')';

K__X_N__X_N = feval(models{g}.covFcn{:}, models{g}.hyp.cov, X_N', []);

%% Mean prediction

% evaluate mean vector
% TODO: tohle asi neni dobre!!! -- nemelo by to byt z X_N ?!
ms = feval(models{g}.meanFcn, models{g}.hyp.mean, X_star');
% noise variance of likGauss
sn2 = exp(2*models{g}.hyp.lik);

% Cholesky factor of covariance with noise
L = chol(K__X_N__X_N/sn2+eye(N) + 0.0001*eye(N));
% Covariance * y_N
K_inv = solve_chol(L,eye(N));

% Predictive mean
Fmu = ms + K__X_star__X_N * K_inv * (1/sn2) * (y_N - ms);

%% Select one of the points
s = 1;
withoutS = [1:(s-1), (s+1):lambda]';
Y_s = X_star(:,s);

X_star_m = X_star(:,withoutS);

% covariances
K__X_star_m__X_N_p = feval(models{g}.covFcn{:}, models{g}.hyp.cov, [X_N X_star(:,s)]', X_star_m')';
K__X_N_p__X_N_p = feval(models{g}.covFcn{:}, models{g}.hyp.cov, [X_N X_star(:,s)]', []);
% evaluate mean vector
% TODO: tohle asi neni dobre!!! -- nemelo by to byt z X_N ?!
ms_m = feval(models{g}.meanFcn, models{g}.hyp.mean, [X_star Y_s]');
% noise variance of likGauss
sn2 = exp(2*models{g}.hyp.lik);
% Cholesky factor of covariance with noise
Lp = chol(K__X_N_p__X_N_p/sn2+eye(N+1) + 0.0001*eye(N+1));
% Covariance * y_N
Kp_inv = solve_chol(Lp,eye(N+1));

% Predictive mean
Fmu_m = ms_m + K__X_star_m__X_N_p * Kp_inv * (1/sn2) * ([y_N  - ms_m);


%% 
A = K__X_star_m__X_N_p * K_inv * (1/sn2);

withoutS = [1:(s-1), (s+1):lambda];
M = NaN(lambda-1,lambda-1);

for i = 1:length(withoutS)
  for j = (i+1):length(withoutS)
    k = withoutS(i);
    l = wihtoutS(j);
    v1 = 
    
  end  
end

