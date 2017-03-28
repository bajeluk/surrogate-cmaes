modelOptions.useShift        = false;
modelOptions.predictionType  = 'sd2';
modelOptions.trainAlgorithm  = 'fmincon';
modelOptions.covFcn          = '{@covMaterniso, 5}';
modelOptions.normalizeY      = true;
modelOptions.hyp.lik         = log(0.01);
modelOptions.hyp.cov         = log([0.5; 2]);
modelOptions.covBounds       = [ [-2;-2], [25;25] ];
modelOptions.likBounds       = log([1e-6, 10]);

% Full factorial design of the following parameters
modelOptions.trainsetType    = { 'nearestToPopulation', 'nearest', 'clustering', 'allPoints' };
modelOptions.trainRange      = { 1.0, 0.999 };
modelOptions.trainsetSizeMax = { '5*dim', '10*dim', '15*dim', '20*dim' };
modelOptions.meanFcn         = { 'meanConst', 'meanLinear' };
modelOptions.covFcn          = { '{@covSEiso}', '{@covSEard}', ...
                            '{@covMaterniso, 5}', '{@covMaterniso, 3}' };
