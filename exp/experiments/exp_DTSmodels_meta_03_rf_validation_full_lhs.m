% Script containing settings of 400 RF models performance testing on
% generation samples from DTS-CMA-ES runs stored in file 
% DTS_meta_005_validation.mat.
% The experiment differs from exp_DTSmodels_meta_03_rf_validation_lhs.m in
% training set selection method here set to 'full', i.e., the model is
% trained on the whole archive.
% RF models differ in decision tree splitting method settings (CART, OC1, 
% PAIR, SECRET, SUPPORT), number of trees in RF (64, 128, 256, 512, 1024),
% number of tree points ({0.25, 0.5, 0.75, 1}*number of available points), 
% and number of tree dimensions ({0.25, 0.5, 0.75, 1}*number of available 
% dimensions). Individual settings are created using full-factorial design.
% However, only 100 settings selected using lhs design will be tested. The
% numbers of selected setting combinations can be found in shell script for
% metacentrum run exp_DTSmodels_meta_03_rf_validation_full_lhs.sh
% Values were taken from pitra2018boosted.

modelOptions.predictionType = 'poi'; % irrelevant - necessary only when opts.testModelOutput = true;
modelOptions.normalizeY = 1;
modelOptions.restartDesign = 'normal';
modelOptions.useShift = 0;
modelOptions.trainsetType = 'full';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

modelOptions.forestType = 'xgboost';
modelOptions.rf_boosting = true;
modelOptions.tree_growFull = false;
modelOptions.tree_predictorFunc = @ConstantModel;
modelOptions.tree_minParentSize = 4;
modelOptions.weak_modelSpec = 'constant';
modelOptions.tree_splitGainFunc = @GradientSplitGain;
modelOptions.tree_maxDepth = 8;
modelOptions.splitGain_modelFunc = @Constant;
modelOptions.split_nQuantize = 10;
modelOptions.split_nRepeats = 1000;
modelOptions.split_maxHyp = '10*dim';

modelOptions.tree_splitFunc = { @AxisSplit, @GaussianSplit, @HillClimbingObliqueSplit, @PairObliqueSplit, @ResidualObliqueSplit };
modelOptions.rf_nTrees = num2cell(2.^[6, 7, 8, 9, 10]);
modelOptions.rf_nFeaturesToSample = { 'ceil(0.25*dim)', 'ceil(0.5*dim)', 'ceil(0.75*dim)', 'dim' };
modelOptions.rf_inBagFraction = num2cell( [0.25, 0.5, 0.75, 1.0] );

opts.snapshotsToTest    = 1:25;

opts.alwaysRetrain      = true;
opts.trySecondModel     = true;
opts.statistics         = { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde', 'rde2', 'rde2models', 'rdeValid', 'rdeValid2', 'rdeM1_M1WReplace', 'rdeM1_M2WReplace', 'rdeM2_M2WReplace', 'mae', 'r2' };
opts.testOrigRatio      = 0.05;
opts.dataset            = 'exp/experiments/dataset/DTS_meta_005_validation.mat';
opts.saveModels         = false;
opts.modelType = 'forest';

% opts.aggFunction        = @(x) nanmedian(reshape(x, [], 1));
% opts.aggSnapshots       = { 1:3, 4:6 };

opts.scratch = [getenv('SCRATCH') '/tmp/' getenv('USER')];
[~, ~] = mkdir(opts.scratch);