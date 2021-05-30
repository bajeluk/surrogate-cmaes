% Script containing settings for calculation of metafeatures on
% generation samples from DTS-CMA-ES runs stored in file
% DTS_meta_005_train.mat.
%
% Unlike in exp_DTSmodels_meta_03_mfts_train, TSS 'knn' with k set
% according to Auger et al. (2013) is calculated. Archive features would be
% identical, thus their computation is not performed.

modelOptions.predictionType = 'poi';
modelOptions.trainAlgorithm = 'fmincon';
modelOptions.normalizeY = 1;
modelOptions.restartDesign = 'normal';
modelOptions.meanFcn = 'meanConst';
modelOptions.useShift = 0;
modelOptions.normalizeY = 1;
modelOptions.trainsetType = 'knn';
modelOptions.knn = 'ceil(min(2 * (obj.dim * (obj.dim + 3)/2 + 1), sqrt(nData * (obj.dim * (obj.dim + 3)/2 + 1))))';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

% no model settings required
opts.modelType = 'random'; % this field is compulsory however irrelevant

opts.snapshotsToTest    = 1:25;
opts.dataset            = 'exp/experiments/dataset/DTS_meta_005_train.mat';
opts.mfts_only          = true;

opts.scratch = [getenv('SCRATCH') '/tmp/' getenv('USER')];
[~, ~] = mkdir(opts.scratch);

% metafeatures settings
opts.mfts_settings.lb = 'min(X)';
opts.mfts_settings.ub = 'max(X)';
opts.mfts_settings.features = {'basic', 'cmaes', 'dispersion', 'ela_distribution', 'ela_levelset', 'ela_metamodel', 'infocontent', 'nearest_better'};
opts.mfts_settings.MetaInput = {'train',     'train', 'traintest',  'traintest'};
opts.mfts_settings.TransData = { 'none',       'cma',      'none',        'cma'};
opts.mfts_settings.useFeat   = {    1:8, [2, 3, 5:8], [1:3, 5, 7], [2, 3, 5, 7]};
opts.mfts_settings.trainOpts = modelOptions;
opts.mfts_settings.warnings = false;
opts.mfts_settings.mixTrans = false; % features are present in two different transformations => do not combine
opts.mfts_settings.Statistics = {'values'};
opts.mfts_settings.StatsOnly = true;
opts.mfts_settings.usePar = false;
opts.mfts_settings.output = 'exp/experiments/exp_DTSmodels_meta_03_mfts_train_knn/metafeatures';
opts.mfts_settings.infocontent.nan_state = true; % consider NaN's as valid value in y
%% use when only metafeatures are required
% opts.mfts_settings.rewrite = true;

%% run
% opts.exp_id = 'exp_DTSmodels_meta_03_mfts_train_knn';
% testModels(modelOptions, opts, 1:24, [2, 3, 5, 10, 20], 11:15, 1:7)
