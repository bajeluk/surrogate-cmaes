modelOptions.useShift        = false;
modelOptions.predictionType  = 'sd2';
modelOptions.normalizeY      = true;
modelOptions.modelType = 'modelPool';
opts.modelType = 'modelPool';

modelOptions.parameterSets = {struct('covFcn', {'{@covSEard}', '{@covSEard}', '{@covSEard}', '{@covSEiso}', '{@covSEard}', '{@covSEard}', '{@covSEiso}', '{@covSEiso}', '{@covSEard}', '{@covMaterniso, 3}', '{@covMaterniso, 3}', '{@covMaterniso, 5}', '{@covMaterniso, 3}', '{@covMaterniso, 5}', '{@covMaterniso, 3}'},'trainsetType', {'nearest', 'nearest', 'nearest', 'nearest', 'allPoints', 'clustering', 'clustering', 'allPoints', 'allPoints', 'nearest', 'allPoints', 'allPoints', 'clustering', 'clustering', 'nearestToPopulation'},'trainRange', {1, 1, 1, 0.999, 1, 1, 1, 0.999, 0.999, 1, 1, 1, 1, 1, 1},'trainsetSizeMax', {'5*dim', '10*dim', '20*dim', '10*dim', '5*dim', '5*dim', '20*dim', '5*dim', '15*dim', '5*dim', '5*dim', '20*dim', '10*dim', '15*dim', '10*dim'},'meanFcn', {'meanLinear', 'meanLinear', 'meanLinear', 'meanConst', 'meanLinear', 'meanLinear', 'meanLinear', 'meanLinear', 'meanLinear', 'meanConst', 'meanConst', 'meanLinear', 'meanConst', 'meanConst', 'meanConst'},'trainAlgorithm', {'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon', 'fmincon'},'hyp', {struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147]), struct('lik', -4.605170, 'cov', [-0.693147, 0.693147])})};

% Full factorial design of the following parameters                          
modelOptions.bestModelSelection = { 'rdeAll', 'rdeOrig' };
modelOptions.historyLength      = { 3, 5, 7 };
modelOptions.minTrainedModelsPercentileForModelChoice = {0.25, 0.5};
modelOptions.maxGenerationShiftForModelChoice = {0, 2};
                          