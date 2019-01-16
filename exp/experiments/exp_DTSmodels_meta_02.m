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

% ADD
modelOptions.hypOptions{end+1}               = struct();
modelOptions.hypOptions{end}.covFcn          = '{@covADD, {[1 max(2, ceil(obj.dim/2))], @covSEisoU}}';
modelOptions.hypOptions{end}.hyp             = struct('lik', log(0.1), 'cov', { [log(0.5) 1 1] } );
modelOptions.hypOptions{end}.hypRestartSigma = 'diag([2 2 * ones(1, obj.dim) 2 2 0.01])';

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
opts.dataset            = 'exp/experiments/dataset/DTS_meta_005.mat';
opts.saveModels         = false;
opts.modelOptionsIndices = [1:6 8:9]; % leaving out ADD kernel for now
opts.modelType = 'gp';

% opts.aggFunction        = @(x) nanmedian(reshape(x, [], 1));
% opts.aggSnapshots       = { 1:3, 4:6 };

opts.scratch = [getenv('SCRATCH') '/tmp/' getenv('USER')];
[~, ~] = mkdir(opts.scratch);

% metafeatures settings
opts.mfts_settings.lb = 'min(X)';
opts.mfts_settings.ub = 'max(X)';
opts.mfts_settings.features = {'basic','cmaes','dispersion','ela_distribution','ela_levelset','ela_metamodel','infocontent','nearest_better'};
opts.mfts_settings.MetaInput = {'archive', 'train', 'train', 'traintest', 'traintest'};
opts.mfts_settings.TransData = {'none',    'none',  'cma',   'none',      'cma'};
opts.mfts_settings.useFeat = {1:8, 2, [1, 3:8], 2, [1, 3:8]};
opts.mfts_settings.trainOpts = modelOptions;
opts.mfts_settings.warnings = false;
opts.mfts_settings.mixTrans = true;
opts.mfts_settings.output = 'exp/experiments/exp_DTSmodels_meta_02/metafeatures';
%% use when only metafeatures are required
% opts.mfts_settings.rewrite = true;
