% Script containing settings of 8 GP models performance testing on
% generation samples from DTS-CMA-ES runs stored in file
% DTS_meta_005_train.mat.
% GP models differ in covariance function settings (LIN, QUAD, SE, MATERN5,
% RQ, NN, SE + QUAD, GIBBS).
% TSS method is 'full', i.e., all archive points are used for training.

modelOptions.predictionType = 'poi';
modelOptions.trainAlgorithm = 'fmincon';
modelOptions.normalizeY = 1;
modelOptions.restartDesign = 'normal';
modelOptions.meanFcn = 'meanConst';
modelOptions.useShift = 0;
modelOptions.normalizeY = 1;
modelOptions.trainsetType = 'full';
modelOptions.nRestarts = 2;
modelOptions.cmaesCheckBounds = false;

modelOptions.hypOptions = {};

% LIN
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covPoly, ''eye'', 1}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.1 1]) });
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 0.01]);

% QUAD
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covPoly, ''eye'', 2}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.1 1]) });
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 0.01]);

% SE
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '@covSEiso';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.5 2]) } );
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 0.01]);

% MATERN5
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covMaterniso, 5}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.5 2]) } );
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 0.01]);

% RQ
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '@covRQiso';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.5 2 0.1]) } );
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 5 0.01]);

% NN (ARCSIN)
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '@covNNone';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.5 2]) } );
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 5 5 0.01]);

% SE + QUAD
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covSum, {@covSEiso, {@covPoly, ''eye'', 2}}}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { log([0.5 2 0.1 1]) } );
modelOptions.hypOptions{end}.hypRestartSigma = diag([2 2 2 2 2 0.01]);

% GIBBS
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covSEvlen, {@meanSum, {@meanLinear, @meanConst}}}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.01), 'cov', { '[zeros(1, obj.dim) 0.5 log([2])]' } );
modelOptions.hypOptions{end}.hypRestartSigma = 'diag([2 10 * ones(1, obj.dim) 10 5 0.01])';

opts.snapshotsToTest    = 1:25;

opts.alwaysRetrain      = true;
opts.trySecondModel     = true;
opts.statistics         = { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde', 'rde2', 'rde2models', 'rdeValid', 'rdeValid2', 'rdeM1_M1WReplace', 'rdeM1_M2WReplace', 'rdeM2_M2WReplace', 'mae', 'r2' };
opts.testOrigRatio      = 0.05;
opts.dataset            = 'exp/experiments/dataset/DTS_meta_005_train.mat';
opts.saveModels         = false;
opts.modelType = 'gp';

opts.scratch = [getenv('SCRATCH') '/tmp/' getenv('USER')];
[~, ~] = mkdir(opts.scratch);