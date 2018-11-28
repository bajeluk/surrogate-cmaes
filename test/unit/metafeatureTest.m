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
  
  % test all NaN values
  y = NaN(nData, 1);
  ft = feature_linear_model(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
  
  % test empty input
  X = [];
  y = [];
  ft = feature_linear_model(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_gcm(X, y);
  % feature_gcm with empty input does not return NaN - is that correct?
%   featFields = fieldnames(ft);
%   for m = 1:numel(featFields)
%     verifyTrue(testCase, isnan(ft.(featFields{m})))
%   end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_cm_convexity(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_cm_gradhomo(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_cm_angle(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_infocontent(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  for m = 1:numel(featFields)
    verifyGreaterThanOrEqual(testCase, ft.(featFields{m}), 0)
    verifyLessThanOrEqual(testCase, ft.(featFields{m}), 1)
  end
  
  % test NaN values
  y = [y(1:20, :); NaN(10, 1)];
  ft = feature_pca(X, y);
  for m = 1:numel(featFields)
    if isempty(strfind(featFields{m}, 'init'))
      verifyTrue(testCase, ~isnan(ft.(featFields{m})))
    else
      verifyTrue(testCase, isnan(ft.(featFields{m})))
    end
  end
  % test all NaN values
  y = NaN(30, 1);
  ft = feature_pca(X, y);
  for m = 1:numel(featFields)
    if isempty(strfind(featFields{m}, 'init'))
      verifyTrue(testCase, ~isnan(ft.(featFields{m})))
    else
      verifyTrue(testCase, isnan(ft.(featFields{m})))
    end
  end
  
  % test empty input
  X = [];
  y = [];
  ft = feature_pca(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_nearest_better(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_dispersion(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test NaN values
  y = [y(1:5); NaN(5, 1)];
  ft = feature_ela_distribution(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, ~isnan(ft.(featFields{m})))
  end
  
  % test full NaN input
  y = NaN(10, 1);
  ft = feature_ela_distribution(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
  
  % test empty input
  X = [];
  y = [];
  ft = feature_ela_distribution(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
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
  
  % test empty input
  X = [];
  y = [];
  ft = feature_ela_levelset(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
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
  
  % test full NaN input
  y = NaN(30, 1);
  ft = feature_ela_metamodel(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
  
  % test empty input
  X = [];
  y = [];
  ft = feature_ela_metamodel(X, y);
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
end

function testBasic(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_basic()));
  % test data
  dim = 10;
  X = rand(30, dim);
  y = randn(30, 1);
  settings.lb = zeros(1, dim);
  settings.ub = ones(1, dim);
  settings.blocks = randi(3, 1, dim);
  % output fields without settings
  featFields = {'dim', 'observations', 'lower_min', 'lower_max', 'upper_min', ...
                'upper_max', 'objective_min', 'objective_max', 'blocks_min', ...
                'blocks_max', 'cells_total', 'cells_filled', 'minimize_fun'};
  returnedFields = fieldnames(feature_basic(X, y));
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  
  % output with settings
  ft = feature_basic(X, y, settings);
  for m = 1:numel(featFields)
    verifyTrue(testCase, ~isnan(ft.(featFields{m})))
  end
  
  % test empty input
  X = [];
  y = [];
  feature_basic(X, y);
end

function testCMA(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, fieldnames(feature_cmaes()));
  % test data
  X = [-2.025294, 1.578529, -2.106864, -3.601003, -2.289839;
       -5, 1.219689, -1.481228, -5, -3.095459; 
       2.882095, 3.838150, -2.825522, -2.757357, 0.691988; 
       -0.100039, 5, -3.633954, -4.906271, -2.461931; 
       0.907036, 0.820051, 2.657243, -1.655711, -0.0418342; 
       -2.299987, 5, -2.629150, -3.800688, -3.776006; 
       -3.163940, -0.676862, 2.277342, -3.303794, -1.588325;
       -1.231081, 2.684261, -2.799170, -5, -3.700963; 
       -2.055466, 1.893202, -3.591501, -4.783686, 0.846982; 
       -2.735506, 3.224880, 0.277076, 0.823197, 3.848332; 
       4.755088, 2.633280, -3.945182, -5, -4.919588; 
       0.483944, 2.535650, 2.210185, -5, -4.094849; 
       -1.033765, 0.394731, -3.273230, -5, -5; 
       -3.212749, -3.161436, -5, -5, -4.136176; 
       -1.648272, 3.299578, -5, 3.644712, -2.113251; 
       -0.490984, -0.204822, -3.174846, -2.164550, -3.019452; 
       2.482224, 0.161941, 1.501581, -2.607810, -1.309724; 
       -3.609463, 3.160187, -1.722410, -5, 1.423379; 
       -0.740302, -3.238282, 0.2152, 2.134075, -5; 
       2.780738, 3.605220, 2.198055, 4.744246, -2.235266; 
       5, 0.0872025, -0.798288, 5, 5; 
       5, -5, 4.671962, 5, -0.410233; 
       3.304363, -5, -3.591884, 5, -5; 
       -4.331972, -1.686711, -2.801912, 4.932170, -5; 
       -0.316491, -4.138878, 4.863647, 1.454858, -5; 
       -2.072590, 2.963153, -1.568964, 4.419250, -5; 
       -1.196841, -4.288049, 1.535959, 4.676978, -1.410081; 
       1.240085, -0.719945, -2.505620, 3.372307, -1.174521; 
       1.133255, -3.854700, -3.809328, 1.797134, -1.159012; 
       1.607137, -2.492724, -2.177072, -2.242048, -1.189113; 
       -1.209478, -3.568769, -0.621612, 0.0587944, 0.110906; 
       -1.568672, 0.183439, 0.327225, 2.303160, -4.337546; 
       1.307934, 1.270794, -3.398229, 3.793819, -2.639051; 
       1.283894, -2.254053, 0.0735176, 3.615344, -2.569520; 
       0.165977, -0.541847, -0.737339, 0.88572, -3.888108; 
       1.142923, 0.310447, 0.109431, 3.740282, -3.156699; 
       0.645641, -1.796558, -0.79374, 5, -4.650979; 
       0.465051, -1.740637, -1.394242, 4.667791, -5; 
       1.274771, -3.576657, -1.718702, 3.904767, -4.602655; 
       1.731575, -2.547548, -0.805452, 5, -5; 
       1.762469, -2.566620, -0.656888, 1.370594, -1.876161; 
       0.0261739, -2.637725, -1.198748, 0.646246, -0.263361; 
       1.891009, -1.273296, -0.859838, 3.510664, -5; 
       -0.807895, 0.344398, -1.773292, 1.228365, -3.651532; 
       -1.171576, -0.463601, -1.278378, 1.461359, 0.0823982; 
       -0.0580385, -1.978201, -1.902091, 0.350758, -1.624527; 
       1.197027, -0.925229, -0.927442, 1.161695, -1.754044; 
       -0.0488689, -1.588799, -1.652942, 1.612937, -1.626606];
  y = [124.769083; 161.439962; 149.072393; 172.711806; 115.044425; 
    161.531515; 128.940939; 149.757730; 159.803659; 152.457874; 
    177.478509; 151.751605; 143.394776; 163.885905; 124.512268; 
    103.796681; 113.581530; 178.849583; 91.102733; 125.224467; 
    172.010674; 160.504993; 126.612686; 119.512006; 125.519717; 
    114.168406; 105.673999; 88.179266; 99.385662; 104.811427; 
    98.727907; 88.585830; 97.127022; 85.248067; 82.406490; 86.636386; 
    93.377103; 93.209055; 94.977030; 98.433293; 84.707563; 89.432692; 
    90.084333; 85.389256; 90.148234; 85.237580; 81.910230; 81.830146];
  BD = [0.233936, -0.135612, -0.2451, 0.168699, -0.159823; 
    -0.00464268, -0.304863, 0.0405761, -0.31218, 0.093645; 
    -0.256917, -0.0929742, -0.266166, 0.11819, -0.0419549; 
    -0.045735, -0.127559, 0.222661, 0.200064, -0.660851; 
    -0.00479693, -0.12268, 0.140304, 0.291703, 0.66289];
  settings.cma_cov = BD*BD;
  settings.cma_evopath_c = [-0.100072; 0.275547; 0.00793752; -1.043735; 0.982928]';
  settings.cma_evopath_s = [0.217826; 0.611778; 0.0882615; -1.217561; 1.058476]';
  settings.cma_generation = 20;
  settings.cma_mean = [0.34444; -1.276915; -0.988959; 1.849136; -2.300293]';
  settings.cma_restart = 1;
  settings.cma_step_size = 1.296304;
  % output fields without settings
  featFields = {'cma_generation', 'cma_step_size', 'cma_restart', ...
                'cma_mean_dist', 'cma_evopath_c_norm', ...
                'cma_evopath_s_norm', ... 'cma_cov_dist', ...
               };
  ft = feature_cmaes(X, y);
  returnedFields = fieldnames(ft);
  verifyTrue(testCase, all(cellfun(@(x) any(strcmp(x, returnedFields)), featFields)))
  % no setting input should cause NaN output
  for m = 1:numel(featFields)
    verifyTrue(testCase, isnan(ft.(featFields{m})))
  end
  % output with settings
  ft = feature_cmaes(X, y, settings);
  for m = 1:numel(featFields)
    verifyTrue(testCase, ~isnan(ft.(featFields{m})))
  end
  
  % test empty input
  X = [];
  y = [];
  feature_cmaes(X, y, settings);
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
  settings.cm_opts.blockType = 'quantile';
  % feature groups
  featGroups =   {'basic', ...
                  'cm_angle', ...
                  'cm_convexity', ...
                  'cm_gradhomo', ...
                  'cmaes', ...
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
  settings.features = {'cm_angle', 'cm_convexity', 'cm_gradhomo'};
  settings.blocks = 3*ones(1, dim);
  settings.cm_gradhomo.blocks = 3*ones(1, dim);
  tic
  [mf, values] = getMetaFeatures(X, y, settings);
  toc
  % two feature groups
  verifyEqual(testCase, numel(struct2cell(mf)), 2)
  % six features
  verifyEqual(testCase, numel(values), 6)
end

function testGetDataMetaFeatures(testCase)
  % empty input should not generate error or warning
  verifyWarningFree(testCase, @getDataMetaFeatures);
  % test data
  testdata = 'exp/experiments/test/DTS_meta';
  outputData = 'exp/experiments/test/DTS_meta_test_fts';
  if ~isdir(testdata)
    warning('Could not finish testGetDataMetaFeatures due to missing test files in %s', ...
            testdata)
    return
  end
  % feature groups
  featGroups =   {'basic', ...
                  'cm_angle', ...
                  'cm_convexity', ...
                  'cm_gradhomo', ...
                  'cmaes', ...
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
  outputFolder = [testdata, '_fts'];
%   getGetMetaFeatures(testdata);
%   verifyTrue(testCase, isdir(outputFolder))
  % test with settings
  settings.lb = 'min(X)'; %'-5*ones(1, dim)';
  settings.ub = 'max(X)'; % 5*ones(1, dim)';
  settings.features = {'basic', 'ela_distribution'}; ... {'cmaes', 'cm_convexity', 'cm_gradhomo'};
  settings.MetaInput = {'traintest'}; % {'archive', 'traintest'};
  settings.trainOpts.evoControlTrainNArchivePoints = '15*dim';
  settings.trainOpts.evoControlTrainRange = 10;
  settings.trainOpts.trainRange = 4;
  settings.trainOpts.trainsetSizeMax = '20*dim';
  settings.trainOpts.trainsetType = 'nearest';
  settings.output = outputData;
  settings.transData = 'cma';
  tic
  getDataMetaFeatures(testdata, settings);
  toc
  verifyTrue(testCase, isdir(outputFolder))
end