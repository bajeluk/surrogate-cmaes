modelOptions.useShift   = 0;
modelOptions.predictionType = 'poi';
modelOptions.trainAlgorithm = 'fmincon';
modelOptions.normalizeY = 1;
modelOptions.meanFcn    = 'meanConst';
modelOptions.hyp.lik    = log(0.01);
modelOptions.hyp.cov    = log([0.5; 2]);
modelOptions.covBounds  = [ [-2;-2], [25;25] ];
modelOptions.likBounds  = log([1e-6, 10]);

modelOptions.trainsetType    = { 'nearestToPopulation', 'nearest', 'clustering', 'recent' };
modelOptions.trainRange      = { 1.5, 2.0, 4.0 };
modelOptions.trainsetSizeMax = { '5*dim', '10*dim', '15*dim', '20*dim' };
modelOptions.covFcn          = { '{@covSEiso}', '{@covMaterniso, 5}', '{@covMaterniso, 3}' };

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
