exp_id = 'exp_doubleEC_MP_18_noisy';
exp_description = 'MP, 10D, noisy, version 13-24 - final_model_settings.m, Copy of DTC_23 - Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, 2-10D, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop, noisy functions';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 10 }, ...
  'functions',          num2cell(101:106), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          num2cell( [1:5, 41:50] ), ...    % default is [1:5, 41:50]
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
structure = struct();

structure(1).covFcn = '{@covMaterniso, 3}';
structure(1).trainsetType = 'allPoints';
structure(1).trainRange = 1;
structure(1).trainsetSizeMax = '5*dim';
structure(1).meanFcn = 'meanConst';
structure(1).trainAlgorithm = 'fmincon';
structure(1).hyp.lik = -4.605170;
structure(1).hyp.cov = [-0.693147; 0.693147];
structure(2).covFcn = '{@covMaterniso, 5}';
structure(2).trainsetType = 'allPoints';
structure(2).trainRange = 1;
structure(2).trainsetSizeMax = '10*dim';
structure(2).meanFcn = 'meanConst';
structure(2).trainAlgorithm = 'fmincon';
structure(2).hyp.lik = -4.605170;
structure(2).hyp.cov = [-0.693147; 0.693147];
structure(3).covFcn = '{@covMaterniso, 5}';
structure(3).trainsetType = 'allPoints';
structure(3).trainRange = 1;
structure(3).trainsetSizeMax = '15*dim';
structure(3).meanFcn = 'meanConst';
structure(3).trainAlgorithm = 'fmincon';
structure(3).hyp.lik = -4.605170;
structure(3).hyp.cov = [-0.693147; 0.693147];
structure(4).covFcn = '{@covMaterniso, 5}';
structure(4).trainsetType = 'allPoints';
structure(4).trainRange = 0.999;
structure(4).trainsetSizeMax = '15*dim';
structure(4).meanFcn = 'meanConst';
structure(4).trainAlgorithm = 'fmincon';
structure(4).hyp.lik = -4.605170;
structure(4).hyp.cov = [-0.693147; 0.693147];
structure(5).covFcn = '{@covMaterniso, 5}';
structure(5).trainsetType = 'nearestToPopulation';
structure(5).trainRange = 1;
structure(5).trainsetSizeMax = '15*dim';
structure(5).meanFcn = 'meanConst';
structure(5).trainAlgorithm = 'fmincon';
structure(5).hyp.lik = -4.605170;
structure(5).hyp.cov = [-0.693147; 0.693147];

modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentileForModelChoice', {0.5},...
    'maxGenerationShiftForModelChoice', {2},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', {  structure  }};



% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
