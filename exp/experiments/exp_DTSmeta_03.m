% Experiment for calculation of metafeature statistics on results from
% fixed experiment exp_doubleEC_28_log_nonadapt using artificial archive
% generation.
%
% Unlike in exp_DTSmeta_02, features are calculated on both transformed and
% non-transformed sets. New set containing archive and testing points is
% added. NaN values are treated in infocontent features. Slightly different
% stats are calculated.

modelOptions.restartDesign = 'normal';
modelOptions.useShift = 0;
modelOptions.normalizeY = 1;
modelOptions.trainsetType = 'nearest';
modelOptions.trainRange = 4.0;
modelOptions.trainsetSizeMax = '20*dim';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

opts.nSampleArchives = 100;
opts.nSnapshotsPerRun = 100;
opts.isForFeatures = true;
opts.isForData = false;
opts.nSmoothGenerations = 2;
opts.smoothDistribution = 'linear';
opts.outputDirname = 'exp_DTSmeta_03';

% metafeatures settings
opts.mfts_settings.lb = 'min(X)';
opts.mfts_settings.ub = 'max(X)';
opts.mfts_settings.features = {'basic', 'cmaes', 'dispersion', 'ela_distribution', 'ela_levelset', 'ela_metamodel', 'infocontent', 'nearest_better'};
opts.mfts_settings.MetaInput = {'archive',   'archive', 'archivetest', 'archivetest', 'train',     'train', 'traintest',  'traintest'};
opts.mfts_settings.TransData = {   'none',       'cma',        'none',         'cma',  'none',       'cma',      'none',        'cma'};
opts.mfts_settings.useFeat   = {      1:8, [2, 3, 5:8],   [1:3, 5, 7],  [2, 3, 5, 7],     1:8, [2, 3, 5:8], [1:3, 5, 7], [2, 3, 5, 7]};
opts.mfts_settings.trainOpts = modelOptions;
opts.mfts_settings.warnings = false;
opts.mfts_settings.mixTrans = false; % features are present in two different transformations => do not combine
opts.mfts_settings.Statistics = {'mean', 'median', 'nan', 'inf', 'var', 'values'};
opts.mfts_settings.StatsOnly = true;
opts.mfts_settings.usePar = false;
opts.mfts_settings.output = 'exp/experiments/exp_DTSmeta_03/sampled_metafeatures';
opts.mfts_settings.infocontent.nan_state = true; % consider NaN's as valid value in y
%% use when only metafeatures are required
% opts.mfts_settings.rewrite = true;

%% run
% modelTestSets('exp_doubleEC_28_log_nonadapt', 1:24, [2, 3, 5, 10, 20], 11:15, opts);