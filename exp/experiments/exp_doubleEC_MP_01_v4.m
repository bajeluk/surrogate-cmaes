exp_id = 'exp_doubleEC_MP_01_v4';
exp_description = 'DTS-CMA-ES with 4 settings of ModelPools: without ARD covariance functions & #3 and #4 also without 5*dim trainsetSizeMax, in 2D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2 }, ...
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
%

% 1) w/o ARD and also w/o 5*dim

% == Best set of settings #1 (id = 1)== 

parameterSets_1a5(1).covFcn = '{@covSEiso}';
parameterSets_1a5(1).trainsetType = 'allPoints';
parameterSets_1a5(1).trainRange = 1;
parameterSets_1a5(1).trainsetSizeMax = '15*dim';
parameterSets_1a5(1).meanFcn = 'meanLinear';
parameterSets_1a5(1).trainAlgorithm = 'fmincon';
parameterSets_1a5(1).hyp.lik = -4.605170;
parameterSets_1a5(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(2).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(2).trainsetType = 'nearest';
parameterSets_1a5(2).trainRange = 1;
parameterSets_1a5(2).trainsetSizeMax = '20*dim';
parameterSets_1a5(2).meanFcn = 'meanConst';
parameterSets_1a5(2).trainAlgorithm = 'fmincon';
parameterSets_1a5(2).hyp.lik = -4.605170;
parameterSets_1a5(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(3).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(3).trainsetType = 'nearest';
parameterSets_1a5(3).trainRange = 1;
parameterSets_1a5(3).trainsetSizeMax = '15*dim';
parameterSets_1a5(3).meanFcn = 'meanConst';
parameterSets_1a5(3).trainAlgorithm = 'fmincon';
parameterSets_1a5(3).hyp.lik = -4.605170;
parameterSets_1a5(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(4).covFcn = '{@covMaterniso, 3}';
parameterSets_1a5(4).trainsetType = 'allPoints';
parameterSets_1a5(4).trainRange = 1;
parameterSets_1a5(4).trainsetSizeMax = '10*dim';
parameterSets_1a5(4).meanFcn = 'meanConst';
parameterSets_1a5(4).trainAlgorithm = 'fmincon';
parameterSets_1a5(4).hyp.lik = -4.605170;
parameterSets_1a5(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(5).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(5).trainsetType = 'allPoints';
parameterSets_1a5(5).trainRange = 1;
parameterSets_1a5(5).trainsetSizeMax = '15*dim';
parameterSets_1a5(5).meanFcn = 'meanLinear';
parameterSets_1a5(5).trainAlgorithm = 'fmincon';
parameterSets_1a5(5).hyp.lik = -4.605170;
parameterSets_1a5(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(6).covFcn = '{@covMaterniso, 3}';
parameterSets_1a5(6).trainsetType = 'clustering';
parameterSets_1a5(6).trainRange = 1;
parameterSets_1a5(6).trainsetSizeMax = '20*dim';
parameterSets_1a5(6).meanFcn = 'meanLinear';
parameterSets_1a5(6).trainAlgorithm = 'fmincon';
parameterSets_1a5(6).hyp.lik = -4.605170;
parameterSets_1a5(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(7).covFcn = '{@covMaterniso, 3}';
parameterSets_1a5(7).trainsetType = 'clustering';
parameterSets_1a5(7).trainRange = 0.999;
parameterSets_1a5(7).trainsetSizeMax = '20*dim';
parameterSets_1a5(7).meanFcn = 'meanConst';
parameterSets_1a5(7).trainAlgorithm = 'fmincon';
parameterSets_1a5(7).hyp.lik = -4.605170;
parameterSets_1a5(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(8).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(8).trainsetType = 'nearestToPopulation';
parameterSets_1a5(8).trainRange = 1;
parameterSets_1a5(8).trainsetSizeMax = '10*dim';
parameterSets_1a5(8).meanFcn = 'meanConst';
parameterSets_1a5(8).trainAlgorithm = 'fmincon';
parameterSets_1a5(8).hyp.lik = -4.605170;
parameterSets_1a5(8).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #2 (id = 8)== 

parameterSets_8a5(1).covFcn = '{@covSEiso}';
parameterSets_8a5(1).trainsetType = 'nearest';
parameterSets_8a5(1).trainRange = 0.999;
parameterSets_8a5(1).trainsetSizeMax = '15*dim';
parameterSets_8a5(1).meanFcn = 'meanLinear';
parameterSets_8a5(1).trainAlgorithm = 'fmincon';
parameterSets_8a5(1).hyp.lik = -4.605170;
parameterSets_8a5(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_8a5(2).covFcn = '{@covMaterniso, 5}';
parameterSets_8a5(2).trainsetType = 'nearest';
parameterSets_8a5(2).trainRange = 1;
parameterSets_8a5(2).trainsetSizeMax = '10*dim';
parameterSets_8a5(2).meanFcn = 'meanLinear';
parameterSets_8a5(2).trainAlgorithm = 'fmincon';
parameterSets_8a5(2).hyp.lik = -4.605170;
parameterSets_8a5(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_8a5(3).covFcn = '{@covMaterniso, 3}';
parameterSets_8a5(3).trainsetType = 'allPoints';
parameterSets_8a5(3).trainRange = 1;
parameterSets_8a5(3).trainsetSizeMax = '10*dim';
parameterSets_8a5(3).meanFcn = 'meanConst';
parameterSets_8a5(3).trainAlgorithm = 'fmincon';
parameterSets_8a5(3).hyp.lik = -4.605170;
parameterSets_8a5(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_8a5(4).covFcn = '{@covMaterniso, 5}';
parameterSets_8a5(4).trainsetType = 'allPoints';
parameterSets_8a5(4).trainRange = 1;
parameterSets_8a5(4).trainsetSizeMax = '10*dim';
parameterSets_8a5(4).meanFcn = 'meanConst';
parameterSets_8a5(4).trainAlgorithm = 'fmincon';
parameterSets_8a5(4).hyp.lik = -4.605170;
parameterSets_8a5(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_8a5(5).covFcn = '{@covMaterniso, 3}';
parameterSets_8a5(5).trainsetType = 'allPoints';
parameterSets_8a5(5).trainRange = 1;
parameterSets_8a5(5).trainsetSizeMax = '15*dim';
parameterSets_8a5(5).meanFcn = 'meanConst';
parameterSets_8a5(5).trainAlgorithm = 'fmincon';
parameterSets_8a5(5).hyp.lik = -4.605170;
parameterSets_8a5(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_8a5(6).covFcn = '{@covMaterniso, 3}';
parameterSets_8a5(6).trainsetType = 'clustering';
parameterSets_8a5(6).trainRange = 0.999;
parameterSets_8a5(6).trainsetSizeMax = '20*dim';
parameterSets_8a5(6).meanFcn = 'meanConst';
parameterSets_8a5(6).trainAlgorithm = 'fmincon';
parameterSets_8a5(6).hyp.lik = -4.605170;
parameterSets_8a5(6).hyp.cov = [-0.693147; 0.693147];

% 2) w/o ARD (5*dim is allowed)

% == Best set of settings #1 (id = 1)== 

parameterSets_1a(1).covFcn = '{@covSEiso}';
parameterSets_1a(1).trainsetType = 'allPoints';
parameterSets_1a(1).trainRange = 1;
parameterSets_1a(1).trainsetSizeMax = '5*dim';
parameterSets_1a(1).meanFcn = 'meanConst';
parameterSets_1a(1).trainAlgorithm = 'fmincon';
parameterSets_1a(1).hyp.lik = -4.605170;
parameterSets_1a(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(2).covFcn = '{@covSEiso}';
parameterSets_1a(2).trainsetType = 'allPoints';
parameterSets_1a(2).trainRange = 1;
parameterSets_1a(2).trainsetSizeMax = '5*dim';
parameterSets_1a(2).meanFcn = 'meanLinear';
parameterSets_1a(2).trainAlgorithm = 'fmincon';
parameterSets_1a(2).hyp.lik = -4.605170;
parameterSets_1a(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(3).covFcn = '{@covSEiso}';
parameterSets_1a(3).trainsetType = 'allPoints';
parameterSets_1a(3).trainRange = 0.999;
parameterSets_1a(3).trainsetSizeMax = '5*dim';
parameterSets_1a(3).meanFcn = 'meanLinear';
parameterSets_1a(3).trainAlgorithm = 'fmincon';
parameterSets_1a(3).hyp.lik = -4.605170;
parameterSets_1a(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(4).covFcn = '{@covMaterniso, 5}';
parameterSets_1a(4).trainsetType = 'allPoints';
parameterSets_1a(4).trainRange = 1;
parameterSets_1a(4).trainsetSizeMax = '5*dim';
parameterSets_1a(4).meanFcn = 'meanConst';
parameterSets_1a(4).trainAlgorithm = 'fmincon';
parameterSets_1a(4).hyp.lik = -4.605170;
parameterSets_1a(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(5).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(5).trainsetType = 'allPoints';
parameterSets_1a(5).trainRange = 1;
parameterSets_1a(5).trainsetSizeMax = '10*dim';
parameterSets_1a(5).meanFcn = 'meanConst';
parameterSets_1a(5).trainAlgorithm = 'fmincon';
parameterSets_1a(5).hyp.lik = -4.605170;
parameterSets_1a(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(6).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(6).trainsetType = 'clustering';
parameterSets_1a(6).trainRange = 1;
parameterSets_1a(6).trainsetSizeMax = '20*dim';
parameterSets_1a(6).meanFcn = 'meanLinear';
parameterSets_1a(6).trainAlgorithm = 'fmincon';
parameterSets_1a(6).hyp.lik = -4.605170;
parameterSets_1a(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(7).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(7).trainsetType = 'clustering';
parameterSets_1a(7).trainRange = 0.999;
parameterSets_1a(7).trainsetSizeMax = '20*dim';
parameterSets_1a(7).meanFcn = 'meanConst';
parameterSets_1a(7).trainAlgorithm = 'fmincon';
parameterSets_1a(7).hyp.lik = -4.605170;
parameterSets_1a(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(8).covFcn = '{@covSEiso}';
parameterSets_1a(8).trainsetType = 'nearestToPopulation';
parameterSets_1a(8).trainRange = 1;
parameterSets_1a(8).trainsetSizeMax = '10*dim';
parameterSets_1a(8).meanFcn = 'meanLinear';
parameterSets_1a(8).trainAlgorithm = 'fmincon';
parameterSets_1a(8).hyp.lik = -4.605170;
parameterSets_1a(8).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(9).covFcn = '{@covMaterniso, 5}';
parameterSets_1a(9).trainsetType = 'nearestToPopulation';
parameterSets_1a(9).trainRange = 0.999;
parameterSets_1a(9).trainsetSizeMax = '5*dim';
parameterSets_1a(9).meanFcn = 'meanConst';
parameterSets_1a(9).trainAlgorithm = 'fmincon';
parameterSets_1a(9).hyp.lik = -4.605170;
parameterSets_1a(9).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #2 (id = 2)== 

parameterSets_2a(1).covFcn = '{@covSEiso}';
parameterSets_2a(1).trainsetType = 'allPoints';
parameterSets_2a(1).trainRange = 1;
parameterSets_2a(1).trainsetSizeMax = '5*dim';
parameterSets_2a(1).meanFcn = 'meanConst';
parameterSets_2a(1).trainAlgorithm = 'fmincon';
parameterSets_2a(1).hyp.lik = -4.605170;
parameterSets_2a(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(2).covFcn = '{@covSEiso}';
parameterSets_2a(2).trainsetType = 'clustering';
parameterSets_2a(2).trainRange = 1;
parameterSets_2a(2).trainsetSizeMax = '20*dim';
parameterSets_2a(2).meanFcn = 'meanConst';
parameterSets_2a(2).trainAlgorithm = 'fmincon';
parameterSets_2a(2).hyp.lik = -4.605170;
parameterSets_2a(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(3).covFcn = '{@covMaterniso, 5}';
parameterSets_2a(3).trainsetType = 'allPoints';
parameterSets_2a(3).trainRange = 1;
parameterSets_2a(3).trainsetSizeMax = '5*dim';
parameterSets_2a(3).meanFcn = 'meanConst';
parameterSets_2a(3).trainAlgorithm = 'fmincon';
parameterSets_2a(3).hyp.lik = -4.605170;
parameterSets_2a(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(4).covFcn = '{@covMaterniso, 3}';
parameterSets_2a(4).trainsetType = 'allPoints';
parameterSets_2a(4).trainRange = 1;
parameterSets_2a(4).trainsetSizeMax = '10*dim';
parameterSets_2a(4).meanFcn = 'meanLinear';
parameterSets_2a(4).trainAlgorithm = 'fmincon';
parameterSets_2a(4).hyp.lik = -4.605170;
parameterSets_2a(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(5).covFcn = '{@covMaterniso, 3}';
parameterSets_2a(5).trainsetType = 'clustering';
parameterSets_2a(5).trainRange = 0.999;
parameterSets_2a(5).trainsetSizeMax = '20*dim';
parameterSets_2a(5).meanFcn = 'meanConst';
parameterSets_2a(5).trainAlgorithm = 'fmincon';
parameterSets_2a(5).hyp.lik = -4.605170;
parameterSets_2a(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(6).covFcn = '{@covMaterniso, 3}';
parameterSets_2a(6).trainsetType = 'nearestToPopulation';
parameterSets_2a(6).trainRange = 1;
parameterSets_2a(6).trainsetSizeMax = '5*dim';
parameterSets_2a(6).meanFcn = 'meanConst';
parameterSets_2a(6).trainAlgorithm = 'fmincon';
parameterSets_2a(6).hyp.lik = -4.605170;
parameterSets_2a(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_2a(7).covFcn = '{@covMaterniso, 5}';
parameterSets_2a(7).trainsetType = 'nearestToPopulation';
parameterSets_2a(7).trainRange = 0.999;
parameterSets_2a(7).trainsetSizeMax = '5*dim';
parameterSets_2a(7).meanFcn = 'meanConst';
parameterSets_2a(7).trainAlgorithm = 'fmincon';
parameterSets_2a(7).hyp.lik = -4.605170;
parameterSets_2a(7).hyp.cov = [-0.693147; 0.693147];


modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', {0.25},...
    'maxGenerationShiftForModelChoice', {2},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_1a5, parameterSets_8a5, parameterSets_1a, parameterSets_2a },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
