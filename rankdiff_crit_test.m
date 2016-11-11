load 'exp/experiments/exp_doubleEC_11_test/bbob_output/exp_doubleEC_11_test_modellog_8_2D_1.mat'
load 'exp/experiments/exp_doubleEC_11_test/exp_doubleEC_11_test_results_8_2D_1.mat'

g = 10;
last_point_idx = cmaes_out{1}{1}.generationStarts(g+1)-1;

orig_120 = logical(cmaes_out{1}{1}.origEvaled(1:last_point_idx));

% X_N = cmaes_out{1}{1}.arxvalids(:,find(orig_120));
% y_N = cmaes_out{1}{1}.fvalues(find(orig_120));
% y_N = y_N';

% TODO: use model's dataset instead!

X_star = cmaes_out{1}{1}.arxvalids(:,cmaes_out{1}{1}.generationStarts(g):last_point_idx);
[f_star, cov_star] = models{10}.predict(X_star');

K__X_star__X_N = feval(models{10}.covFcn{:}, models{10}.hyp.cov, models{10}.dataset.X, X_star')';

K__X_N__X_N = feval(models{10}.covFcn{:}, models{10}.hyp.cov, models{10}.dataset.X, []);

[D N] = size(X_N);

% evaluate mean vector
ms = feval(models{10}.meanFcn, models{10}.hyp.mean, X_star');
% noise variance of likGauss
sn2 = exp(2*models{10}.hyp.lik);
% sqrt of noise precision vector
post_sW = ones(N,1)/sqrt(sn2);

% Cholesky factor of covariance with noise
L = chol(K__X_N__X_N/sn2+eye(N) + 0.0001*eye(N));

alpha = solve_chol(L,y_N'-ms)/sn2;

% number of alphas (usually 1; more in case of sampling)
% TODO: is it equal to N?!
N2 = size(alpha,2);
% conditional mean fs|f
Fmu = repmat(ms,1,N2) + K__X_star__X_N * alpha;
% predictive means
fmu(id) = sum(Fmu,2)/N2;


