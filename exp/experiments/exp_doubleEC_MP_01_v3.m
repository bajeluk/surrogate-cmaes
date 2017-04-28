exp_id = 'exp_doubleEC_MP_01_v3';
exp_description = 'DTS-CMA-ES with 3 shorter settings of ModelPools w/o ARD in 2D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

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
% == Best set of settings #1 (id = 14)== 
parameterSet_14(1).covFcn = '{@covSEiso}';
parameterSet_14(1).trainsetType = 'clustering';
parameterSet_14(1).trainRange = 0.999;
parameterSet_14(1).trainsetSizeMax = '20*dim';
parameterSet_14(1).meanFcn = 'meanLinear';
parameterSet_14(1).trainAlgorithm = 'fmincon';
parameterSet_14(1).hyp.lik = -4.605170;
parameterSet_14(1).hyp.cov = [-0.693147; 0.693147];
parameterSet_14(2).covFcn = '{@covMaterniso, 5}';
parameterSet_14(2).trainsetType = 'allPoints';
parameterSet_14(2).trainRange = 1;
parameterSet_14(2).trainsetSizeMax = '5*dim';
parameterSet_14(2).meanFcn = 'meanConst';
parameterSet_14(2).trainAlgorithm = 'fmincon';
parameterSet_14(2).hyp.lik = -4.605170;
parameterSet_14(2).hyp.cov = [-0.693147; 0.693147];
parameterSet_14(3).covFcn = '{@covMaterniso, 3}';
parameterSet_14(3).trainsetType = 'nearestToPopulation';
parameterSet_14(3).trainRange = 1;
parameterSet_14(3).trainsetSizeMax = '15*dim';
parameterSet_14(3).meanFcn = 'meanConst';
parameterSet_14(3).trainAlgorithm = 'fmincon';
parameterSet_14(3).hyp.lik = -4.605170;
parameterSet_14(3).hyp.cov = [-0.693147; 0.693147];
parameterSet_14(4).covFcn = '{@covMaterniso, 5}';
parameterSet_14(4).trainsetType = 'nearestToPopulation';
parameterSet_14(4).trainRange = 1;
parameterSet_14(4).trainsetSizeMax = '15*dim';
parameterSet_14(4).meanFcn = 'meanConst';
parameterSet_14(4).trainAlgorithm = 'fmincon';
parameterSet_14(4).hyp.lik = -4.605170;
parameterSet_14(4).hyp.cov = [-0.693147; 0.693147];
parameterSet_14(5).covFcn = '{@covMaterniso, 5}';
parameterSet_14(5).trainsetType = 'nearestToPopulation';
parameterSet_14(5).trainRange = 0.999;
parameterSet_14(5).trainsetSizeMax = '5*dim';
parameterSet_14(5).meanFcn = 'meanLinear';
parameterSet_14(5).trainAlgorithm = 'fmincon';
parameterSet_14(5).hyp.lik = -4.605170;
parameterSet_14(5).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #2 (id = 2)== 
parameterSet_2(1).trainRange = 1;
parameterSet_2(1).trainsetSizeMax = '5*dim';
parameterSet_2(1).meanFcn = 'meanLinear';
parameterSet_2(1).trainAlgorithm = 'fmincon';
parameterSet_2(1).hyp.lik = -4.605170;
parameterSet_2(1).hyp.cov = [-0.693147; 0.693147];
parameterSet_2(2).trainRange = 1;
parameterSet_2(2).trainsetSizeMax = '5*dim';
parameterSet_2(2).meanFcn = 'meanLinear';
parameterSet_2(2).trainAlgorithm = 'fmincon';
parameterSet_2(2).hyp.lik = -4.605170;
parameterSet_2(2).hyp.cov = [-0.693147; 0.693147];
parameterSet_2(3).trainRange = 1;
parameterSet_2(3).trainsetSizeMax = '5*dim';
parameterSet_2(3).meanFcn = 'meanConst';
parameterSet_2(3).trainAlgorithm = 'fmincon';
parameterSet_2(3).hyp.lik = -4.605170;
parameterSet_2(3).hyp.cov = [-0.693147; 0.693147];
parameterSet_2(4).covFcn = '{@covMaterniso, 3}';
parameterSet_2(4).trainsetType = 'allPoints';
parameterSet_2(4).trainRange = 1;
parameterSet_2(4).trainsetSizeMax = '10*dim';
parameterSet_2(4).meanFcn = 'meanLinear';
parameterSet_2(4).trainAlgorithm = 'fmincon';
parameterSet_2(4).hyp.lik = -4.605170;
parameterSet_2(4).hyp.cov = [-0.693147; 0.693147];
parameterSet_2(5).covFcn = '{@covMaterniso, 3}';
parameterSet_2(5).trainsetType = 'clustering';
parameterSet_2(5).trainRange = 0.999;
parameterSet_2(5).trainsetSizeMax = '20*dim';
parameterSet_2(5).meanFcn = 'meanConst';
parameterSet_2(5).trainAlgorithm = 'fmincon';
parameterSet_2(5).hyp.lik = -4.605170;
parameterSet_2(5).hyp.cov = [-0.693147; 0.693147];
parameterSet_2(6).trainRange = 0.999;
parameterSet_2(6).trainsetSizeMax = '5*dim';
parameterSet_2(6).meanFcn = 'meanConst';
parameterSet_2(6).trainAlgorithm = 'fmincon';
parameterSet_2(6).hyp.lik = -4.605170;
parameterSet_2(6).hyp.cov = [-0.693147; 0.693147];

% == Best set of settings #3 (id = 1)== 
parameterSet_1(1).trainRange = 0.999;
parameterSet_1(1).trainsetSizeMax = '20*dim';
parameterSet_1(1).meanFcn = 'meanConst';
parameterSet_1(1).trainAlgorithm = 'fmincon';
parameterSet_1(1).hyp.lik = -4.605170;
parameterSet_1(1).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(2).covFcn = '{@covMaterniso, 5}';
parameterSet_1(2).trainsetType = 'nearest';
parameterSet_1(2).trainRange = 1;
parameterSet_1(2).trainsetSizeMax = '5*dim';
parameterSet_1(2).meanFcn = 'meanLinear';
parameterSet_1(2).trainAlgorithm = 'fmincon';
parameterSet_1(2).hyp.lik = -4.605170;
parameterSet_1(2).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(3).covFcn = '{@covMaterniso, 5}';
parameterSet_1(3).trainsetType = 'allPoints';
parameterSet_1(3).trainRange = 1;
parameterSet_1(3).trainsetSizeMax = '5*dim';
parameterSet_1(3).meanFcn = 'meanConst';
parameterSet_1(3).trainAlgorithm = 'fmincon';
parameterSet_1(3).hyp.lik = -4.605170;
parameterSet_1(3).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(4).covFcn = '{@covMaterniso, 3}';
parameterSet_1(4).trainsetType = 'allPoints';
parameterSet_1(4).trainRange = 1;
parameterSet_1(4).trainsetSizeMax = '10*dim';
parameterSet_1(4).meanFcn = 'meanLinear';
parameterSet_1(4).trainAlgorithm = 'fmincon';
parameterSet_1(4).hyp.lik = -4.605170;
parameterSet_1(4).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(5).covFcn = '{@covMaterniso, 3}';
parameterSet_1(5).trainsetType = 'clustering';
parameterSet_1(5).trainRange = 0.999;
parameterSet_1(5).trainsetSizeMax = '20*dim';
parameterSet_1(5).meanFcn = 'meanConst';
parameterSet_1(5).trainAlgorithm = 'fmincon';
parameterSet_1(5).hyp.lik = -4.605170;
parameterSet_1(5).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(6).trainRange = 0.999;
parameterSet_1(6).trainsetSizeMax = '20*dim';
parameterSet_1(6).meanFcn = 'meanLinear';
parameterSet_1(6).trainAlgorithm = 'fmincon';
parameterSet_1(6).hyp.lik = -4.605170;
parameterSet_1(6).hyp.cov = [-0.693147; 0.693147];
parameterSet_1(7).covFcn = '{@covMaterniso, 5}';
parameterSet_1(7).trainsetType = 'nearestToPopulation';
parameterSet_1(7).trainRange = 0.999;
parameterSet_1(7).trainsetSizeMax = '5*dim';
parameterSet_1(7).meanFcn = 'meanConst';
parameterSet_1(7).trainAlgorithm = 'fmincon';
parameterSet_1(7).hyp.lik = -4.605170;
parameterSet_1(7).hyp.cov = [-0.693147; 0.693147];


modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', {0.25},...
    'maxGenerationShiftForModelChoice', {2},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSet_1, parameterSet_2, parameterSet_14 },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
