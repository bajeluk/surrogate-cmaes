%Original experiment: exp_doubleEC_26_1model
exp_id = 'marek_lmmEC_lmm_validation';
exp_description = 'lmm-CMA-ES original, in [2,3,5,10,20]D, on instances [1:5, 31:40], PreSampleSize=0.75, 1pop';

% BBOB/COCO framework settings

bbobParams = { ...
  'dimensions',         { 2, 3, 5, 10, 20 }, ...
  'functions',          num2cell(1:24), ...      % all functions: num2cell(1:24)
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          { [1:5, 31:40] }, ...    % default is [1:5, 41:50]
  'maxfunevals',        { '250 * dim' }, ...
  'resume',             { true }, ...
};

% Surrogate manager parameters

surrogateParams = { ...
  'evoControl',         { 'lmm' }, ...    % 'none', 'individual', 'generation', 'restricted'
  'observers',          { {'DTScreenStatistics', 'DTFileStatistics'} },... % logging observers
  'modelType',          { 'lmm' }, ...               % 'gp', 'rf', 'bbob'
  'updaterType',        { 'rankDiff' }, ...         % OrigRatioUpdater
  'evoControlMaxDoubleTrainIterations', { 1 }, ...
  'evoControlPreSampleSize',       { 0.75 }, ...       % {0.25, 0.5, 0.75}, will be multip. by lambda
  'evoControlOrigPointsRoundFcn',  { 'ceil' }, ...  % 'ceil', 'getProbNumber'
  'evoControlTrainRange',          { 100000 }, ...      % will be multip. by sigma
  'evoControlSampleRange',         { 1 }, ...       % will be multip. by sigma
  'evoControlOrigGenerations',     { [] }, ...
  'evoControlModelGenerations',    { [] }, ...
  'evoControlValidatePoints',      { [] }, ...
  'evoControlRestrictedParam',     { 0.05 }, ...
  'evoControlUseInject',           { false }, ...
};

% Model parameters

modelParams = { ...
  'covFcn',             { '{@covMaterniso, 5}' }, ...
  'hyp',                { struct('lik', log(0.01), 'cov', log([0.5; 2])) }, ...
  'meanFcn',            { 'meanConst' }, ...
  'trainAlgorithm',     { 'fmincon' }, ...
  'predictionType',     { 'fvalues' }, ...
  'useShift',           { false }, ...
  'normalizeY',         { true }, ...
  'trainsetType'        { 'nearest' }, ...
  'trainRange',         { 100000 }, ...
  'trainsetSizeMax'     { '20*dim' }, ...l
};


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(4+floor(3*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/hanusm12/public';
