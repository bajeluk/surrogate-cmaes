function tests = metafeatureTest
% unit test for metafeature functions
  tests = functiontests(localfunctions);
end

function testLinearModel(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_linear_model()));
  % test data
  dim = 20;
  nData = 50*dim;
  X = rand(nData, dim);
  y = randn(nData, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  settings.blocks = [4, 4, 3*ones(1, dim-2)];
  % output fields with settings
  featFields = {'lm_avg_length_reg', 'lm_avg_length_norm', 'lm_length_mean', ...
                'lm_length_std', 'lm_corr_reg', 'lm_corr_norm', 'lm_ratio_mean', ...
                'lm_ratio_std', 'lm_std_radio_reg', 'lm_std_radio_norm', ...
                'lm_std_mean_reg', 'lm_std_mean_norm'};
  ft = feature_linear_model(X, y, settings);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  for f = 1 : numel(returnedFields)
    verifySize(testCase, ft.(returnedFields{f}), [1 1])
  end
end

function testGCM(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_gcm()));
  % test data
  dim = 5;
  nData = 50*dim;
  X = rand(nData, dim);
  y = randn(nData, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  % settings.blocks = [2*ones(1, dim)];
  settings.blocks = [4, 4, 3*ones(1, dim-2)];
  % output fields without settings
  ft = feature_gcm(X, y, settings);
  returnedFields = fieldnames(ft);
  nRetFields = numel(returnedFields);
  % number of returned fields can be divided by 23 and be at least 23
  verifyTrue(testCase, mod(nRetFields, 23) == 0 && nRetFields >= 23)
  % output range
  for m = 1:nRetFields
    % verifyGreaterThanOrEqual(testCase, ft.(returnedFields{m}), 0)
    % if ~strfind(ft.(returnedFields{m}), 'attractor')
    %   verifyLessThanOrEqual(testCase, ft.(returnedFields{m}), 1)
    % end
  end
end

function testCMConvexity(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_cm_convexity()));
  % test data
  dim = 20;
  nData = 50*dim;
  X = rand(nData, dim);
  y = randn(nData, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  settings.blocks = [4, 4, 3*ones(1, dim-2)];
  % output fields with settings
  featFields = {'concave_soft', 'concave_hard', 'convex_soft', 'convex_hard'};
  ft = feature_cm_convexity(X, y, settings);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testCMGradHomo(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_cm_gradhomo()));
  % test data
  dim = 10;
  nData = 300;
  X = rand(nData, dim);
  y = randn(nData, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  settings.blocks = 2*ones(1, dim);
  % output fields with settings
  featFields = {'grad_mean', 'grad_std'};
  ft = feature_cm_gradhomo(X, y, settings);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testCMAngle(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_cm_angle()));
  % test data
  dim = 20;
  nPoints = 50*dim;
  X = rand(nPoints, dim);
  y = randn(nPoints, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  settings.blocks = randi(3, 1, dim);
  % output fields with settings
  featFields = {'dist_ctr2best_mean', 'dist_ctr2best_std', ...
                'dist_ctr2worst_mean', 'dist_ctr2worst_std', ...
                'angle_mean', 'angle_std', ...
                'y_best2worst_mean', 'y_best2worst_std'};
  ft = feature_cm_angle(X, y, settings);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testInfocontent(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_infocontent()));
  % test data
  X = rand(30, 10);
  y = randn(30, 1);
  % output fields without settings
  featFields = {'h_max', 'eps_s', 'eps_max', 'eps_ratio', 'm0'};
  ft = feature_infocontent(X, y);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  % test random sorting
  settings.sorting = 'random';
  ft = feature_infocontent(X, y, settings);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testPCA(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_pca()));
  % test data
  X = rand(30, 10);
  y = randn(30, 1);
  % output fields without settings
  featFields = {'pca_cov_x', 'pca_corr_x', 'pca_cov_init', 'pca_corr_init', ...
                'pca_pc1_cov_x', 'pca_pc1_corr_x', 'pca_pc1_cov_init', ...
                'pca_pc1_corr_init'};
  ft = feature_pca(X, y);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  % test percentage values
  for m = 1:4
    verifyGreaterThanOrEqual(testCase, ft.(featFields{m}), 0)
    verifyLessThanOrEqual(testCase, ft.(featFields{m}), 1)
  end
end

function testNearestBetterClustering(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_nearest_better()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % test settings
  settings.distance = 'minkowski';
  % output fields with settings
  featFields = {'nb_std_ratio', 'nb_mean_ratio', 'nb_cor', 'dist_ratio', ...
                'nb_fitness_cor'};
  returnedFields = fieldnames(feature_nearest_better(X, y, settings));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testDispersion(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_dispersion()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % test settings
  settings.distance = 'minkowski';
  settings.quantiles = [0.1, 0.25];
  % output fields with settings
  featFields = {'ratio_mean_10', 'ratio_median_10', 'diff_mean_10', 'diff_median_10', ...
                'ratio_mean_25', 'ratio_median_25', 'diff_mean_25', 'diff_median_25'};
  returnedFields = fieldnames(feature_dispersion(X, y, settings));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testELADistribution(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_distribution()));
  % test data
  X = rand(10, 3);
  y = randn(10, 1);
  % output fields
  featFields = {'skewness', 'kurtosis', 'number_of_peaks'};
  returnedFields = fieldnames(feature_ela_distribution(X, y));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testELALevelset(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_levelset()));
  % test data
  nData = 100;
  X = rand(nData, 3);
  y = randn(nData, 1);
  qnt = [10, 25, 50];
  % output fields
  featFields = {'mmce_lda_10', 'mmce_lda_25', 'mmce_lda_50', ...
                'mmce_qda_10', 'mmce_qda_25', 'mmce_qda_50', ... 
                'lda_qda_10',  'lda_qda_25',  'lda_qda_50', ...
                'mmce_mda_10', 'mmce_mda_25', 'mmce_mda_50', ...
                'lda_mda_10',  'lda_mda_25',  'lda_mda_50', ...
                'qda_mda_10',  'qda_mda_25',  'qda_mda_50'};
  ft = feature_ela_levelset(X, y);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  % test mmce values
  for m = 1:numel(returnedFields)
    if ~isempty(strfind(returnedFields{m}, 'mmce'))
      verifyGreaterThanOrEqual(testCase, ft.(returnedFields{m}), 0)
      verifyLessThanOrEqual(testCase, ft.(returnedFields{m}), 1)
    end
  end
  % test combination values
  for q = qnt
    verifyEqual(testCase, ...
                ft.(sprintf('lda_qda_%d', q)), ...
                ft.(sprintf('mmce_lda_%d', q)) / ft.(sprintf('mmce_qda_%d', q)))
  end
end

function testELAMetamodel(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_ela_metamodel()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % output fields
  featFields = {'lin_simple_adj_r2', 'lin_simple_intercept', 'lin_simple_coef_min', ...
                'lin_simple_coef_max', 'lin_simple_coef_max_by_min', ...
                'lin_w_interact_adj_r2', 'quad_simple_adj_r2', ...
                'quad_simple_cond', 'quad_w_interact_adj_r2'};
  returnedFields = fieldnames(feature_ela_metamodel(X, y)); 
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testBasic(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_basic()));
  % test data
  X = rand(30, 3);
  y = randn(30, 1);
  % output fields without settings
  featFields = {'dim', 'observations', 'lower_min', 'lower_max', 'upper_min', ...
                'upper_max', 'objective_min', 'objective_max', 'blocks_min', ...
                'blocks_max', 'cells_total', 'cells_filled', 'minimize_fun'};
  returnedFields = fieldnames(feature_basic(X, y));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
end

function testGetMetaFeatures(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(getMetaFeatures()));
  % test data
  dim = 10;
  nData = 50*dim;
  X = rand(nData, dim);
  y = randn(nData, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  % feature groups
  featGroups =   {'basic', ...
                  'cm_angle', ...
                  'cm_convexity', ...
                  'cm_gradhomo', ...
                  'dispersion', ...
                  'ela_distribution', ...
                  'ela_levelset', ...
                  'ela_metamodel', ...
                  'gcm', ...
                  'infocontent', ...
                  'linear_model', ...
                  'nearest_better', ...
                  'pca' ...
                 };
  % test no settings
  mf = getMetaFeatures(X, y);
  returnedFields = fieldnames(mf);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featGroups)))
  % test with settings
  settings.features = {'cm_convexity', 'cm_gradhomo'};
  settings.blocks = 2*ones(1, dim);
  settings.cm_gradhomo.blocks = 3*ones(1, dim);
  tic
  [mf, values] = getMetaFeatures(X, y, settings);
  toc
  % two feature groups
  verifyEqual(testCase, numel(struct2cell(mf)), 2)
  % six features
  verifyEqual(testCase, numel(values), 6)
end