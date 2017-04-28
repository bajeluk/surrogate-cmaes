exp_id = 'exp_doubleEC_MP_02_v2';
exp_description = 'DTS-CMA-ES with 3 different settings of ModelPools in 5D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 5 }, ...
  'functions',          num2cell([2, 4, 6, 8, 10, 12, 16, 18, 22, 23]), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          num2cell([1:5, 41:50]), ...    % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { true }, ...
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'doubletrained' }, ...    % 'none', 'individual', 'generation', 'restricted'
  'observers',          { {'DTScreenStatistics', 'DTFileStatistics'} },... % logging observers
  'modelType',          { 'modelPool' }, ...               % 'gp', 'rf', 'bbob'
  'updaterType',        { 'rankDiff' }, ...         % OrigRatioUpdater
  'evoControlMaxDoubleTrainIterations', { 1 }, ...
  'evoControlPreSampleSize',       { 0.75 }, ...       % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlOrigPointsRoundFcn',  { 'ceil' }, ...  % 'ceil', 'getProbNumber'
  'evoControlIndividualExtension', { [] }, ...      % will be multip. by lambda
  'evoControlBestFromExtension',   { [] }, ...      % ratio of expanded popul.
  'evoControlTrainRange',          { 10 }, ...      % will be multip. by sigma
  'evoControlTrainNArchivePoints', { '15*dim' },... % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used
  'evoControlSampleRange',         { 1 }, ...       % will be multip. by sigma
  'evoControlOrigGenerations',     { [] }, ...
  'evoControlModelGenerations',    { [] }, ...
  'evoControlValidatePoints',      { [] }, ...
  'evoControlRestrictedParam',     { 0.05 }, ...
};

% Model parameters

% == Best set of settings #1 (id = 13)== 

parameterSets_13(1).covFcn = '{@covSEard}';
parameterSets_13(1).trainsetType = 'nearest';
parameterSets_13(1).trainRange = 1;
parameterSets_13(1).trainsetSizeMax = '15*dim';
parameterSets_13(1).meanFcn = 'meanConst';
parameterSets_13(1).trainAlgorithm = 'fmincon';
parameterSets_13(1).hyp.lik = -4.605170;
parameterSets_13(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(2).covFcn = '{@covSEiso}';
parameterSets_13(2).trainsetType = 'nearest';
parameterSets_13(2).trainRange = 1;
parameterSets_13(2).trainsetSizeMax = '20*dim';
parameterSets_13(2).meanFcn = 'meanLinear';
parameterSets_13(2).trainAlgorithm = 'fmincon';
parameterSets_13(2).hyp.lik = -4.605170;
parameterSets_13(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(3).covFcn = '{@covSEard}';
parameterSets_13(3).trainsetType = 'allPoints';
parameterSets_13(3).trainRange = 1;
parameterSets_13(3).trainsetSizeMax = '5*dim';
parameterSets_13(3).meanFcn = 'meanConst';
parameterSets_13(3).trainAlgorithm = 'fmincon';
parameterSets_13(3).hyp.lik = -4.605170;
parameterSets_13(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(4).covFcn = '{@covMaterniso, 5}';
parameterSets_13(4).trainsetType = 'nearest';
parameterSets_13(4).trainRange = 1;
parameterSets_13(4).trainsetSizeMax = '5*dim';
parameterSets_13(4).meanFcn = 'meanConst';
parameterSets_13(4).trainAlgorithm = 'fmincon';
parameterSets_13(4).hyp.lik = -4.605170;
parameterSets_13(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(5).covFcn = '{@covSEiso}';
parameterSets_13(5).trainsetType = 'clustering';
parameterSets_13(5).trainRange = 0.999;
parameterSets_13(5).trainsetSizeMax = '10*dim';
parameterSets_13(5).meanFcn = 'meanConst';
parameterSets_13(5).trainAlgorithm = 'fmincon';
parameterSets_13(5).hyp.lik = -4.605170;
parameterSets_13(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(6).covFcn = '{@covMaterniso, 5}';
parameterSets_13(6).trainsetType = 'nearest';
parameterSets_13(6).trainRange = 1;
parameterSets_13(6).trainsetSizeMax = '15*dim';
parameterSets_13(6).meanFcn = 'meanLinear';
parameterSets_13(6).trainAlgorithm = 'fmincon';
parameterSets_13(6).hyp.lik = -4.605170;
parameterSets_13(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(7).covFcn = '{@covMaterniso, 3}';
parameterSets_13(7).trainsetType = 'allPoints';
parameterSets_13(7).trainRange = 1;
parameterSets_13(7).trainsetSizeMax = '10*dim';
parameterSets_13(7).meanFcn = 'meanConst';
parameterSets_13(7).trainAlgorithm = 'fmincon';
parameterSets_13(7).hyp.lik = -4.605170;
parameterSets_13(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(8).covFcn = '{@covMaterniso, 5}';
parameterSets_13(8).trainsetType = 'allPoints';
parameterSets_13(8).trainRange = 1;
parameterSets_13(8).trainsetSizeMax = '10*dim';
parameterSets_13(8).meanFcn = 'meanConst';
parameterSets_13(8).trainAlgorithm = 'fmincon';
parameterSets_13(8).hyp.lik = -4.605170;
parameterSets_13(8).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(9).covFcn = '{@covMaterniso, 5}';
parameterSets_13(9).trainsetType = 'allPoints';
parameterSets_13(9).trainRange = 1;
parameterSets_13(9).trainsetSizeMax = '15*dim';
parameterSets_13(9).meanFcn = 'meanConst';
parameterSets_13(9).trainAlgorithm = 'fmincon';
parameterSets_13(9).hyp.lik = -4.605170;
parameterSets_13(9).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(10).covFcn = '{@covMaterniso, 5}';
parameterSets_13(10).trainsetType = 'clustering';
parameterSets_13(10).trainRange = 0.999;
parameterSets_13(10).trainsetSizeMax = '20*dim';
parameterSets_13(10).meanFcn = 'meanConst';
parameterSets_13(10).trainAlgorithm = 'fmincon';
parameterSets_13(10).hyp.lik = -4.605170;
parameterSets_13(10).hyp.cov = [-0.693147; 0.693147];
parameterSets_13(11).covFcn = '{@covMaterniso, 3}';
parameterSets_13(11).trainsetType = 'nearestToPopulation';
parameterSets_13(11).trainRange = 1;
parameterSets_13(11).trainsetSizeMax = '15*dim';
parameterSets_13(11).meanFcn = 'meanConst';
parameterSets_13(11).trainAlgorithm = 'fmincon';
parameterSets_13(11).hyp.lik = -4.605170;
parameterSets_13(11).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #2 (id = 12)== 

parameterSets_12(1).covFcn = '{@covSEard}';
parameterSets_12(1).trainsetType = 'nearest';
parameterSets_12(1).trainRange = 1;
parameterSets_12(1).trainsetSizeMax = '5*dim';
parameterSets_12(1).meanFcn = 'meanConst';
parameterSets_12(1).trainAlgorithm = 'fmincon';
parameterSets_12(1).hyp.lik = -4.605170;
parameterSets_12(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(2).covFcn = '{@covSEiso}';
parameterSets_12(2).trainsetType = 'nearest';
parameterSets_12(2).trainRange = 1;
parameterSets_12(2).trainsetSizeMax = '20*dim';
parameterSets_12(2).meanFcn = 'meanLinear';
parameterSets_12(2).trainAlgorithm = 'fmincon';
parameterSets_12(2).hyp.lik = -4.605170;
parameterSets_12(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(3).covFcn = '{@covSEard}';
parameterSets_12(3).trainsetType = 'allPoints';
parameterSets_12(3).trainRange = 1;
parameterSets_12(3).trainsetSizeMax = '5*dim';
parameterSets_12(3).meanFcn = 'meanConst';
parameterSets_12(3).trainAlgorithm = 'fmincon';
parameterSets_12(3).hyp.lik = -4.605170;
parameterSets_12(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(4).covFcn = '{@covMaterniso, 3}';
parameterSets_12(4).trainsetType = 'allPoints';
parameterSets_12(4).trainRange = 1;
parameterSets_12(4).trainsetSizeMax = '10*dim';
parameterSets_12(4).meanFcn = 'meanConst';
parameterSets_12(4).trainAlgorithm = 'fmincon';
parameterSets_12(4).hyp.lik = -4.605170;
parameterSets_12(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(5).covFcn = '{@covMaterniso, 3}';
parameterSets_12(5).trainsetType = 'nearest';
parameterSets_12(5).trainRange = 0.999;
parameterSets_12(5).trainsetSizeMax = '15*dim';
parameterSets_12(5).meanFcn = 'meanConst';
parameterSets_12(5).trainAlgorithm = 'fmincon';
parameterSets_12(5).hyp.lik = -4.605170;
parameterSets_12(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(6).covFcn = '{@covMaterniso, 5}';
parameterSets_12(6).trainsetType = 'allPoints';
parameterSets_12(6).trainRange = 1;
parameterSets_12(6).trainsetSizeMax = '15*dim';
parameterSets_12(6).meanFcn = 'meanLinear';
parameterSets_12(6).trainAlgorithm = 'fmincon';
parameterSets_12(6).hyp.lik = -4.605170;
parameterSets_12(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_12(7).covFcn = '{@covMaterniso, 5}';
parameterSets_12(7).trainsetType = 'clustering';
parameterSets_12(7).trainRange = 0.999;
parameterSets_12(7).trainsetSizeMax = '20*dim';
parameterSets_12(7).meanFcn = 'meanConst';
parameterSets_12(7).trainAlgorithm = 'fmincon';
parameterSets_12(7).hyp.lik = -4.605170;
parameterSets_12(7).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #3 (id = 15)== 

parameterSets_15(1).covFcn = '{@covSEiso}';
parameterSets_15(1).trainsetType = 'nearest';
parameterSets_15(1).trainRange = 1;
parameterSets_15(1).trainsetSizeMax = '20*dim';
parameterSets_15(1).meanFcn = 'meanLinear';
parameterSets_15(1).trainAlgorithm = 'fmincon';
parameterSets_15(1).hyp.lik = -4.605170;
parameterSets_15(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(2).covFcn = '{@covSEiso}';
parameterSets_15(2).trainsetType = 'allPoints';
parameterSets_15(2).trainRange = 0.999;
parameterSets_15(2).trainsetSizeMax = '20*dim';
parameterSets_15(2).meanFcn = 'meanConst';
parameterSets_15(2).trainAlgorithm = 'fmincon';
parameterSets_15(2).hyp.lik = -4.605170;
parameterSets_15(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(3).covFcn = '{@covMaterniso, 3}';
parameterSets_15(3).trainsetType = 'allPoints';
parameterSets_15(3).trainRange = 1;
parameterSets_15(3).trainsetSizeMax = '10*dim';
parameterSets_15(3).meanFcn = 'meanConst';
parameterSets_15(3).trainAlgorithm = 'fmincon';
parameterSets_15(3).hyp.lik = -4.605170;
parameterSets_15(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(4).covFcn = '{@covMaterniso, 3}';
parameterSets_15(4).trainsetType = 'nearest';
parameterSets_15(4).trainRange = 0.999;
parameterSets_15(4).trainsetSizeMax = '15*dim';
parameterSets_15(4).meanFcn = 'meanConst';
parameterSets_15(4).trainAlgorithm = 'fmincon';
parameterSets_15(4).hyp.lik = -4.605170;
parameterSets_15(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(5).covFcn = '{@covMaterniso, 5}';
parameterSets_15(5).trainsetType = 'allPoints';
parameterSets_15(5).trainRange = 1;
parameterSets_15(5).trainsetSizeMax = '5*dim';
parameterSets_15(5).meanFcn = 'meanLinear';
parameterSets_15(5).trainAlgorithm = 'fmincon';
parameterSets_15(5).hyp.lik = -4.605170;
parameterSets_15(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(6).covFcn = '{@covMaterniso, 5}';
parameterSets_15(6).trainsetType = 'allPoints';
parameterSets_15(6).trainRange = 1;
parameterSets_15(6).trainsetSizeMax = '15*dim';
parameterSets_15(6).meanFcn = 'meanLinear';
parameterSets_15(6).trainAlgorithm = 'fmincon';
parameterSets_15(6).hyp.lik = -4.605170;
parameterSets_15(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(7).covFcn = '{@covMaterniso, 3}';
parameterSets_15(7).trainsetType = 'clustering';
parameterSets_15(7).trainRange = 0.999;
parameterSets_15(7).trainsetSizeMax = '10*dim';
parameterSets_15(7).meanFcn = 'meanConst';
parameterSets_15(7).trainAlgorithm = 'fmincon';
parameterSets_15(7).hyp.lik = -4.605170;
parameterSets_15(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_15(8).covFcn = '{@covMaterniso, 5}';
parameterSets_15(8).trainsetType = 'clustering';
parameterSets_15(8).trainRange = 0.999;
parameterSets_15(8).trainsetSizeMax = '20*dim';
parameterSets_15(8).meanFcn = 'meanConst';
parameterSets_15(8).trainAlgorithm = 'fmincon';
parameterSets_15(8).hyp.lik = -4.605170;
parameterSets_15(8).hyp.cov = [-0.693147; 0.693147];

modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', {0.5},...
    'maxGenerationShiftForModelChoice', {0},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_12, parameterSets_13, parameterSets_15 },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
