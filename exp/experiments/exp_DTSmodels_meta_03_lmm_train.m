% Script containing setting of lmm model performance testing on generation
% samples from DTS-CMA-ES runs stored in file DTS_meta_005_train.mat.
% TSS nearest with settings from bajer2018gaussian is utilized.

modelOptions.predictionType = 'fvalues';
modelOptions.normalizeY = 1;
modelOptions.restartDesign = 'normal';
modelOptions.useShift = 0;
modelOptions.trainsetType = 'nearest';
modelOptions.trainRange = 4.0;
modelOptions.trainsetSizeMax = '20*dim';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

% lmm specific options
modelOptions.lmmKnn = Inf; % using all points returned using TSS nearest

opts.snapshotsToTest    = 1:25;

opts.alwaysRetrain      = true;
opts.trySecondModel     = true;
opts.statistics         = { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde', 'rde2', 'rde2models', 'rdeValid', 'rdeValid2', 'rdeM1_M1WReplace', 'rdeM1_M2WReplace', 'rdeM2_M2WReplace', 'mae', 'r2' };
opts.testOrigRatio      = 0.05;
opts.dataset            = 'exp/experiments/dataset/DTS_meta_005_train.mat';
opts.saveModels         = false;
opts.modelType          = 'lmm';

opts.scratch = [getenv('SCRATCH') '/tmp/' getenv('USER')];
[~, ~] = mkdir(opts.scratch);