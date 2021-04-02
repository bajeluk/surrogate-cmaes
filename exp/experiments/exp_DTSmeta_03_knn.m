% Experiment for calculation of metafeature statistics on results from
% fixed experiment exp_doubleEC_28_log_nonadapt using artificial archive
% generation.
%
% Unlike in exp_DTSmeta_03, TSS 'knn' with k set according to Auger et al.
% (2013) is calculated. Archive features would be identical, thus their
% computation is not performed.

modelOptions.restartDesign = 'normal';
modelOptions.useShift = 0;
modelOptions.normalizeY = 1;
modelOptions.trainsetType = 'knn';
modelOptions.knn = 'ceil(min(2 * (obj.dim * (obj.dim + 3)/2 + 1), sqrt(nData * (obj.dim * (obj.dim + 3)/2 + 1))))';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

opts.nSampleArchives = 100;
opts.nSnapshotsPerRun = 100;
opts.isForFeatures = true;
opts.isForData = false;
opts.nSmoothGenerations = 2;
opts.smoothDistribution = 'linear';
opts.outputDirname = 'exp_DTSmeta_03_knn';

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
opts.mfts_settings.Statistics = {'mean', 'median', 'nan', 'inf', 'var', 'values'};
opts.mfts_settings.StatsOnly = true;
opts.mfts_settings.usePar = false;
opts.mfts_settings.output = fullfile('exp', 'experiments', opts.outputDirname, 'sampled_metafeatures');
opts.mfts_settings.infocontent.nan_state = true; % consider NaN's as valid value in y
%% use when only metafeatures are required
% opts.mfts_settings.rewrite = true;

%% run
% modelTestSets('exp_doubleEC_28_log_nonadapt', 1:24, [2, 3, 5, 10, 20], 11:15, opts);
