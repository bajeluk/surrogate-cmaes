% exp_GPtest_03
%
% very short EXAMPLE experiment examining four different GP models

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
modelOptions.trainsetType    = { 'nearestToPopulation' };
modelOptions.trainRange      = { 0.999 };
modelOptions.trainsetSizeMax = { '10*dim' };
modelOptions.meanFcn         = { 'meanConst' };
modelOptions.covFcn          = { '{@covSEiso}', '{@covSEard}', ...
                            '{@covMaterniso, 5}', '{@covMaterniso, 3}' };

