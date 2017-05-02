exp_id = 'exp_doubleEC_MP_02_v5';
exp_description = 'DTS-CMA-ES with 1 settings of ModelPools: w/o ARD, and w/o meanLinear, in 5D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

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
%

% w/o ARD and w/o meanLinear (total # = 96)

% rdeRankingOpts.maxRank = {20, 25};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };
% rdeRankingOpts.includeARD = { false };
% rdeRankingOpts.include5dim = { true };
% rdeRankingOpts.includeMeanLinear = { false };

% == Best set of settings #1 (id = 10)== 
%    best according to the max_covrd_RDE
%
%   id    maxRank    colStat    minTrainedPerc    colSetCover                           f_weight                           includeARD    include5dim    includeMeanLinear    nSettings    avg_nCovered    avg_RDE    avg_covrd_RDE    max_covrd_RDE    avg_rank    avg_covrd_rank    avg_n_covrd_funs    n_ARD    n_Linear
%   __    _______    _______    ______________    ___________    ______________________________________________________    __________    ___________    _________________    _________    ____________    _______    _____________    _____________    ________    ______________    ________________    _____    ________
%
%   10    25         2           0.7              1              @(nCover,modelErrors)(nCover.^(4))./sum(modelErrors,2)    false         true           false                5              17.8          0.41625    0.41402          0.43903          39.354      10.884              14.2              0        0       
%
%
%   settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                              covrd_funs                               covrd_snp           covFcn           trainsetType    trainRange    trainsetSizeMax      meanFcn  
%   _________    __________    _______    _________    ________    __________    _______________________________________________________________    _________    ____________________    ____________    __________    _______________    ___________
%
%    6           12            0.42546    0.43903        50.5      9.5833        '2   3   7   8   9  10  11  14  19  21  22'                        10     2     '{@covSEiso}'           'nearest'       0.999         '5*dim'            'meanConst'
%   25           18            0.41121    0.42107      39.208      9.7778        '1   2   3   4   6   8   9  10  11  13  14  16  22  23  24'        11     7     '{@covMaterniso, 5}'    'nearest'           1         '10*dim'           'meanConst'
%   37           25            0.39977    0.38379      27.438           9        '1   2   3   4   5   6   8   9  10  11  12  13  14  19  20  24'    12    13     '{@covMaterniso, 3}'    'allPoints'         1         '10*dim'           'meanConst'
%   47           17            0.41358    0.40751      38.333      11.588        '1   3   4   5  10  11  12  13  14  15  18  20  22  23  24'        13     4     '{@covMaterniso, 3}'    'nearest'       0.999         '15*dim'           'meanConst'
%   68           17            0.43124    0.41871      41.292      14.471        '3   4   6   7  12  15  16  17  18  19  21  22  23  24'             6    11     '{@covMaterniso, 5}'    'clustering'    0.999         '10*dim'           'meanConst'

parameterSets_10woAL(1).covFcn = '{@covSEiso}';
parameterSets_10woAL(1).trainsetType = 'nearest';
parameterSets_10woAL(1).trainRange = 0.999;
parameterSets_10woAL(1).trainsetSizeMax = '5*dim';
parameterSets_10woAL(1).meanFcn = 'meanConst';
parameterSets_10woAL(1).trainAlgorithm = 'fmincon';
parameterSets_10woAL(1).hyp.lik = -4.605170;
parameterSets_10woAL(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(2).covFcn = '{@covMaterniso, 5}';
parameterSets_10woAL(2).trainsetType = 'nearest';
parameterSets_10woAL(2).trainRange = 1;
parameterSets_10woAL(2).trainsetSizeMax = '10*dim';
parameterSets_10woAL(2).meanFcn = 'meanConst';
parameterSets_10woAL(2).trainAlgorithm = 'fmincon';
parameterSets_10woAL(2).hyp.lik = -4.605170;
parameterSets_10woAL(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(3).covFcn = '{@covMaterniso, 3}';
parameterSets_10woAL(3).trainsetType = 'allPoints';
parameterSets_10woAL(3).trainRange = 1;
parameterSets_10woAL(3).trainsetSizeMax = '10*dim';
parameterSets_10woAL(3).meanFcn = 'meanConst';
parameterSets_10woAL(3).trainAlgorithm = 'fmincon';
parameterSets_10woAL(3).hyp.lik = -4.605170;
parameterSets_10woAL(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(4).covFcn = '{@covMaterniso, 3}';
parameterSets_10woAL(4).trainsetType = 'nearest';
parameterSets_10woAL(4).trainRange = 0.999;
parameterSets_10woAL(4).trainsetSizeMax = '15*dim';
parameterSets_10woAL(4).meanFcn = 'meanConst';
parameterSets_10woAL(4).trainAlgorithm = 'fmincon';
parameterSets_10woAL(4).hyp.lik = -4.605170;
parameterSets_10woAL(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(5).covFcn = '{@covMaterniso, 5}';
parameterSets_10woAL(5).trainsetType = 'clustering';
parameterSets_10woAL(5).trainRange = 0.999;
parameterSets_10woAL(5).trainsetSizeMax = '10*dim';
parameterSets_10woAL(5).meanFcn = 'meanConst';
parameterSets_10woAL(5).trainAlgorithm = 'fmincon';
parameterSets_10woAL(5).hyp.lik = -4.605170;
parameterSets_10woAL(5).hyp.cov = [-0.693147; 0.693147];

modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', { 0.5 },...
    'maxGenerationShiftForModelChoice', { 0 },...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_10woAL },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
