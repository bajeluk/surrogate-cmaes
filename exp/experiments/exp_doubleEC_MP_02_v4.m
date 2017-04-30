exp_id = 'exp_doubleEC_MP_02_v4';
exp_description = 'DTS-CMA-ES with 2 settings of ModelPools: without ARD covariance functions & also without 5*dim trainsetSizeMax, in 5D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

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
% 1) w/o ARD and also w/o 5*dim

% == Best set of settings #1 (id = 1)== 
% 
%    settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                                    covrd_funs                                     covrd_snp           covFcn           trainsetType    trainRange    trainsetSizeMax      meanFcn   
%    _________    __________    _______    _________    ________    __________    ___________________________________________________________________________    _________    ____________________    ____________    __________    _______________    ____________
% 
%      3           4                NaN          0      87.396       19.75        '4   7  14  23'                                                                 3     1     '{@covSEiso}'           'nearest'           1         '15*dim'           'meanConst' 
%     28           6            0.46328    0.44262      94.646          24        '3   4   7  19  20  21'                                                         1     5     '{@covSEiso}'           'allPoints'     0.999         '10*dim'           'meanLinear'
%     49          29                NaN          0      34.083      13.759        '1   2   3   4   5   6   8   9  10  11  12  13  14  15  17  18  19  20  24'    15    14     '{@covMaterniso, 3}'    'allPoints'         1         '10*dim'           'meanConst' 
%     58          23                NaN          0      50.125      19.391        '1   2   4   8   9  10  11  12  13  14  16  17  18  19  20  22  24'            11    12     '{@covMaterniso, 5}'    'allPoints'         1         '15*dim'           'meanConst' 
%     74           3                NaN          0      100.67      16.667        '5   6  21'                                                                     0     3     '{@covMaterniso, 5}'    'clustering'        1         '10*dim'           'meanConst' 
%    100          20                NaN          0      48.938       13.65        '3   4   6   7  10  11  12  15  16  17  18  19  20  21  22  23  24'             7    13     '{@covMaterniso, 5}'    'clustering'    0.999         '20*dim'           'meanConst' 

parameterSets_1a5(1).covFcn = '{@covSEiso}';
parameterSets_1a5(1).trainsetType = 'nearest';
parameterSets_1a5(1).trainRange = 1;
parameterSets_1a5(1).trainsetSizeMax = '15*dim';
parameterSets_1a5(1).meanFcn = 'meanConst';
parameterSets_1a5(1).trainAlgorithm = 'fmincon';
parameterSets_1a5(1).hyp.lik = -4.605170;
parameterSets_1a5(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(2).covFcn = '{@covSEiso}';
parameterSets_1a5(2).trainsetType = 'allPoints';
parameterSets_1a5(2).trainRange = 0.999;
parameterSets_1a5(2).trainsetSizeMax = '10*dim';
parameterSets_1a5(2).meanFcn = 'meanLinear';
parameterSets_1a5(2).trainAlgorithm = 'fmincon';
parameterSets_1a5(2).hyp.lik = -4.605170;
parameterSets_1a5(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(3).covFcn = '{@covMaterniso, 3}';
parameterSets_1a5(3).trainsetType = 'allPoints';
parameterSets_1a5(3).trainRange = 1;
parameterSets_1a5(3).trainsetSizeMax = '10*dim';
parameterSets_1a5(3).meanFcn = 'meanConst';
parameterSets_1a5(3).trainAlgorithm = 'fmincon';
parameterSets_1a5(3).hyp.lik = -4.605170;
parameterSets_1a5(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(4).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(4).trainsetType = 'allPoints';
parameterSets_1a5(4).trainRange = 1;
parameterSets_1a5(4).trainsetSizeMax = '15*dim';
parameterSets_1a5(4).meanFcn = 'meanConst';
parameterSets_1a5(4).trainAlgorithm = 'fmincon';
parameterSets_1a5(4).hyp.lik = -4.605170;
parameterSets_1a5(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(5).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(5).trainsetType = 'clustering';
parameterSets_1a5(5).trainRange = 1;
parameterSets_1a5(5).trainsetSizeMax = '10*dim';
parameterSets_1a5(5).meanFcn = 'meanConst';
parameterSets_1a5(5).trainAlgorithm = 'fmincon';
parameterSets_1a5(5).hyp.lik = -4.605170;
parameterSets_1a5(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a5(6).covFcn = '{@covMaterniso, 5}';
parameterSets_1a5(6).trainsetType = 'clustering';
parameterSets_1a5(6).trainRange = 0.999;
parameterSets_1a5(6).trainsetSizeMax = '20*dim';
parameterSets_1a5(6).meanFcn = 'meanConst';
parameterSets_1a5(6).trainAlgorithm = 'fmincon';
parameterSets_1a5(6).hyp.lik = -4.605170;
parameterSets_1a5(6).hyp.cov = [-0.693147; 0.693147];


% 2) w/o ARD (5*dim is allowed)

% == Best set of settings #1 (id = 1)== 
%
%    settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                            covrd_funs                             covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn   
%    _________    __________    _______    _________    ________    __________    ___________________________________________________________    _________    ____________________    _____________________    __________    _______________    ____________
%
%      4           3                NaN          0      115.83          16        '4   7  23'                                                     2     1     '{@covSEiso}'           'nearest'                    1         '15*dim'           'meanConst' 
%      9          16                NaN          0      75.875       16.75        '2   3   5   8   9  10  11  13  14  15  17  19  22  24'         9     7     '{@covSEiso}'           'allPoints'                  1         '5*dim'            'meanConst' 
%     65          17            0.44054      0.421      64.458      14.059        '1   5   8  10  11  12  13  14  17  18  19  20  21  24'        10     7     '{@covMaterniso, 3}'    'allPoints'                  1         '5*dim'            'meanConst' 
%     69          22                NaN          0      43.646      11.318        '1   2   3   4   5   6   8   9  10  11  12  13  14  19  20'    11    11     '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst' 
%     80          15            0.48612    0.48666      74.438      20.067        '1   2   4   6   8   9  10  11  12  14  15  18  22  24'        10     5     '{@covMaterniso, 5}'    'nearest'                0.999         '15*dim'           'meanConst' 
%    107           1                NaN          0      162.79           2        '5'                                                             0     1     '{@covMaterniso, 3}'    'clustering'                 1         '5*dim'            'meanLinear'
%    134          16            0.41806    0.41665      66.583      14.562        '3   6   7  12  15  16  17  18  19  20  21  22  24'             4    12     '{@covMaterniso, 5}'    'clustering'             0.999         '20*dim'           'meanConst' 
%    161          20            0.43889    0.42359      56.521        17.9        '1   2   5   8  11  12  13  14  16  17  18  19  20  22  23'    11     9     '{@covMaterniso, 3}'    'nearestToPopulation'        1         '5*dim'            'meanConst' 

parameterSets_1a(1).covFcn = '{@covSEiso}';
parameterSets_1a(1).trainsetType = 'nearest';
parameterSets_1a(1).trainRange = 1;
parameterSets_1a(1).trainsetSizeMax = '15*dim';
parameterSets_1a(1).meanFcn = 'meanConst';
parameterSets_1a(1).trainAlgorithm = 'fmincon';
parameterSets_1a(1).hyp.lik = -4.605170;
parameterSets_1a(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(2).covFcn = '{@covSEiso}';
parameterSets_1a(2).trainsetType = 'allPoints';
parameterSets_1a(2).trainRange = 1;
parameterSets_1a(2).trainsetSizeMax = '5*dim';
parameterSets_1a(2).meanFcn = 'meanConst';
parameterSets_1a(2).trainAlgorithm = 'fmincon';
parameterSets_1a(2).hyp.lik = -4.605170;
parameterSets_1a(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(3).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(3).trainsetType = 'allPoints';
parameterSets_1a(3).trainRange = 1;
parameterSets_1a(3).trainsetSizeMax = '5*dim';
parameterSets_1a(3).meanFcn = 'meanConst';
parameterSets_1a(3).trainAlgorithm = 'fmincon';
parameterSets_1a(3).hyp.lik = -4.605170;
parameterSets_1a(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(4).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(4).trainsetType = 'allPoints';
parameterSets_1a(4).trainRange = 1;
parameterSets_1a(4).trainsetSizeMax = '10*dim';
parameterSets_1a(4).meanFcn = 'meanConst';
parameterSets_1a(4).trainAlgorithm = 'fmincon';
parameterSets_1a(4).hyp.lik = -4.605170;
parameterSets_1a(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(5).covFcn = '{@covMaterniso, 5}';
parameterSets_1a(5).trainsetType = 'nearest';
parameterSets_1a(5).trainRange = 0.999;
parameterSets_1a(5).trainsetSizeMax = '15*dim';
parameterSets_1a(5).meanFcn = 'meanConst';
parameterSets_1a(5).trainAlgorithm = 'fmincon';
parameterSets_1a(5).hyp.lik = -4.605170;
parameterSets_1a(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(6).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(6).trainsetType = 'clustering';
parameterSets_1a(6).trainRange = 1;
parameterSets_1a(6).trainsetSizeMax = '5*dim';
parameterSets_1a(6).meanFcn = 'meanLinear';
parameterSets_1a(6).trainAlgorithm = 'fmincon';
parameterSets_1a(6).hyp.lik = -4.605170;
parameterSets_1a(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(7).covFcn = '{@covMaterniso, 5}';
parameterSets_1a(7).trainsetType = 'clustering';
parameterSets_1a(7).trainRange = 0.999;
parameterSets_1a(7).trainsetSizeMax = '20*dim';
parameterSets_1a(7).meanFcn = 'meanConst';
parameterSets_1a(7).trainAlgorithm = 'fmincon';
parameterSets_1a(7).hyp.lik = -4.605170;
parameterSets_1a(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_1a(8).covFcn = '{@covMaterniso, 3}';
parameterSets_1a(8).trainsetType = 'nearestToPopulation';
parameterSets_1a(8).trainRange = 1;
parameterSets_1a(8).trainsetSizeMax = '5*dim';
parameterSets_1a(8).meanFcn = 'meanConst';
parameterSets_1a(8).trainAlgorithm = 'fmincon';
parameterSets_1a(8).hyp.lik = -4.605170;
parameterSets_1a(8).hyp.cov = [-0.693147; 0.693147];


modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', {0.25},...
    'maxGenerationShiftForModelChoice', {2},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_1a5, parameterSets_1a },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';

