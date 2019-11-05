% experiment for calculation of metafeature statistics on results from 
% experiment exp_doubleEC_28_log_nonadapt using artificial archive
% generation

modelOptions.predictionType = 'poi';
modelOptions.trainAlgorithm = 'fmincon';
modelOptions.normalizeY = 1;
modelOptions.restartDesign = 'normal';
modelOptions.meanFcn = 'meanConst';
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
opts.outputDirname = 'exp_DTSmeta_02';

% metafeatures settings
opts.mfts_settings.lb = 'min(X)';
opts.mfts_settings.ub = 'max(X)';
opts.mfts_settings.features = {'basic','cmaes','dispersion','ela_distribution','ela_levelset','ela_metamodel','infocontent','nearest_better'};
opts.mfts_settings.MetaInput = {'archive', 'train', 'train', 'traintest', 'traintest'};
opts.mfts_settings.TransData = {'none',    'none',  'cma',   'none',      'cma'};
opts.mfts_settings.useFeat = {1:8, 2, [1, 3:8], 2, [1, 3:5, 7:8]};
opts.mfts_settings.trainOpts = modelOptions;
opts.mfts_settings.warnings = false;
opts.mfts_settings.mixTrans = true;
opts.mfts_settings.Statistics = {'mean', 'var', 'hist'};
opts.mfts_settings.usePar = false;
opts.mfts_settings.output = 'exp/experiments/exp_DTSmeta_02/sampled_metafeatures';
%% use when only metafeatures are required
% opts.mfts_settings.rewrite = true;

%% run
% modelTestSets('exp_doubleEC_28_log_nonadapt', 1:24, [2, 3, 5, 10, 20], 11:15, opts);
