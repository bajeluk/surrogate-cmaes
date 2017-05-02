exp_id = 'exp_doubleEC_MP_02_v6';
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

% rdeRankingOpts.maxRank = {25, 35};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };
% rdeRankingOpts.includeARD = { false };
% rdeRankingOpts.include5dim = { true };
% rdeRankingOpts.includeMeanLinear = { false };

% == Best set of settings #1 (id = 10)== 
%    best according to the avg_covrd_RDE
%
%    id    maxRank    colStat    minTrainedPerc    colSetCover                           f_weight                           includeARD    include5dim    includeMeanLinear    nSettings    avg_nCovered    avg_RDE    avg_covrd_RDE    max_covrd_RDE    avg_rank    avg_covrd_rank    avg_n_covrd_funs    n_ARD    n_Linear
%    __    _______    _______    ______________    ___________    ______________________________________________________    __________    ___________    _________________    _________    ____________    _______    _____________    _____________    ________    ______________    ________________    _____    ________
%
%    10    35         2           0.7              1              @(nCover,modelErrors)(nCover.^(4))./sum(modelErrors,2)    false         true           false                5              25.6          0.41118    0.39802          0.41967          35.879      16.515              17.8              0        0       
%     1    25         2           0.7              1              @(nCover,modelErrors)(nCover.^(3))./sum(modelErrors,2)    false         true           false                7            15.429          0.41456    0.40681          0.46273          43.426      11.765            12.429              0        0       
%    12    35         2          0.85              1              @(nCover,modelErrors)(nCover.^(4))./sum(modelErrors,2)    false         true           false                5              21.2          0.41483    0.40705          0.42296          40.492      16.678              16.4              0        0       
%
% == Best set of settings #1 (id = 10)== 
%
%    settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                                    covrd_funs                                     covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn  
%    _________    __________    _______    _________    ________    __________    ___________________________________________________________________________    _________    ____________________    _____________________    __________    _______________    ___________
% 
%    33           27            0.40564    0.39398      37.417      18.407        '1   2   3   4   5   8  10  11  12  13  14  17  18  19  20  21  23  24'        14    13     '{@covMaterniso, 3}'    'allPoints'                  1         '5*dim'            'meanConst'
%    37           30            0.39977    0.38882      27.438      12.733        '1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  19  20  21  24'    14    16     '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst'
%    46           23            0.41463    0.39624      40.292          18        '1   2   3   4   5   7   8   9  10  11  12  13  14  17  19  20  22  24'        10    13     '{@covMaterniso, 5}'    'allPoints'                  1         '15*dim'           'meanConst'
%    70           21            0.42973    0.41967      40.771      15.619        '3   4   6   7  10  12  15  16  17  18  19  20  21  22  23  24'                 7    14     '{@covMaterniso, 5}'    'clustering'             0.999         '20*dim'           'meanConst'
%    81           27            0.40612     0.3914      33.479      17.815        '1   2   5   6   8   9  10  11  12  13  14  16  17  18  19  20  22  23'        13    14     '{@covMaterniso, 3}'    'nearestToPopulation'        1         '5*dim'            'meanConst'

parameterSets_10woAL(1).covFcn = '{@covMaterniso, 3}';
parameterSets_10woAL(1).trainsetType = 'allPoints';
parameterSets_10woAL(1).trainRange = 1;
parameterSets_10woAL(1).trainsetSizeMax = '5*dim';
parameterSets_10woAL(1).meanFcn = 'meanConst';
parameterSets_10woAL(1).trainAlgorithm = 'fmincon';
parameterSets_10woAL(1).hyp.lik = -4.605170;
parameterSets_10woAL(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(2).covFcn = '{@covMaterniso, 3}';
parameterSets_10woAL(2).trainsetType = 'allPoints';
parameterSets_10woAL(2).trainRange = 1;
parameterSets_10woAL(2).trainsetSizeMax = '10*dim';
parameterSets_10woAL(2).meanFcn = 'meanConst';
parameterSets_10woAL(2).trainAlgorithm = 'fmincon';
parameterSets_10woAL(2).hyp.lik = -4.605170;
parameterSets_10woAL(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(3).covFcn = '{@covMaterniso, 5}';
parameterSets_10woAL(3).trainsetType = 'allPoints';
parameterSets_10woAL(3).trainRange = 1;
parameterSets_10woAL(3).trainsetSizeMax = '15*dim';
parameterSets_10woAL(3).meanFcn = 'meanConst';
parameterSets_10woAL(3).trainAlgorithm = 'fmincon';
parameterSets_10woAL(3).hyp.lik = -4.605170;
parameterSets_10woAL(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(4).covFcn = '{@covMaterniso, 5}';
parameterSets_10woAL(4).trainsetType = 'clustering';
parameterSets_10woAL(4).trainRange = 0.999;
parameterSets_10woAL(4).trainsetSizeMax = '20*dim';
parameterSets_10woAL(4).meanFcn = 'meanConst';
parameterSets_10woAL(4).trainAlgorithm = 'fmincon';
parameterSets_10woAL(4).hyp.lik = -4.605170;
parameterSets_10woAL(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_10woAL(5).covFcn = '{@covMaterniso, 3}';
parameterSets_10woAL(5).trainsetType = 'nearestToPopulation';
parameterSets_10woAL(5).trainRange = 1;
parameterSets_10woAL(5).trainsetSizeMax = '5*dim';
parameterSets_10woAL(5).meanFcn = 'meanConst';
parameterSets_10woAL(5).trainAlgorithm = 'fmincon';
parameterSets_10woAL(5).hyp.lik = -4.605170;
parameterSets_10woAL(5).hyp.cov = [-0.693147; 0.693147];


modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', { 0.5 },...
    'maxGenerationShiftForModelChoice', { 2 },...
    'predictionType',     { 'sd2', 'ei' }, ...
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
