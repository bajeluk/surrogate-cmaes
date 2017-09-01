modelOptions.useShift   = 0;
modelOptions.predictionType = 'poi';
modelOptions.trainAlgorithm = 'fmincon';
modelOptions.normalizeY = 1;

modelOptions.covFcn     = {'{@covMaterniso, 5}', 'covSEiso'};
modelOptions.meanFcn    = {'meanConst', 'meanZero'}
modelOptions.covBounds  = {[ [-2;-2], [25;25] ], [ [log(1);-2], [log(1);25] ]};

modelOptions.hyp.lik    = log(0.01);
% modelOptions.hyp.cov    = {log([0.5; 2]);
modelOptions.likBounds  = log([1e-6, 10]);
modelOptions.trainsetType    = 'nearest';
modelOptions.trainRange      = 4.0;
modelOptions.trainsetSizeMax = '20*dim';

opts.snapshotsToTest    = [ 5, 6, 7, 18, 19, 20 ];

opts.alwaysRetrain      = true;
opts.trySecondModel     = true;
opts.statistics         = { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde', 'rde2', 'rde2models', 'rdeValid', 'rdeValid2', 'rdeM1_M1WReplace', 'rdeM1_M2WReplace', 'rdeM2_M2WReplace' };
opts.testOrigRatio      = 0.05;
opts.dataset            = 'DTS_005_25_models';
opts.saveModels         = true;

opts.aggFunction        = @(x) nanmedian(reshape(x, [], 1));
opts.aggSnapshots       = { 1:3, 4:6 };

opts.scratch = [getenv('SCRATCH') '/tmp/bajeluk'];
[~, ~] = mkdir(opts.scratch);
