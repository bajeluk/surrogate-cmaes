exp_id = 'exp_doubleEC_MP_01_v5';
exp_description = 'DTS-CMA-ES with 1 settings of ModelPools: w/o ARD, w/o 5*dim and w/o meanLinear, in 2D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

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

% 1) w/o ARD, w/o 5*dim and w/o meanLinear

% == Best set of settings #1 (id = 3)== 
%
% rdeRankingOpts.maxRank = {25, 35};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };
%
%     settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                                covrd_funs                                 covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn  
%     _________    __________    _______    _________    ________    __________    ___________________________________________________________________    _________    ____________________    _____________________    __________    _______________    ___________
% 
%     19           11            0.44959    0.42817        43.5      14.455        '1   5   6   7  10  12  14  16  24'                                     9     2     '{@covMaterniso, 5}'    'nearest'                    1         '20*dim'           'meanConst'
%     25           26            0.42295    0.40691      24.292      10.846        '1   2   4   5   8   9  10  11  13  14  15  16  18  19  20  23  24'    15    11     '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst'
%     26           23            0.42078    0.40481      24.292      11.261        '1   2   4   5   6   8   9  10  11  13  15  17  18  19  20  23'        13    10     '{@covMaterniso, 5}'    'allPoints'                  1         '10*dim'           'meanConst'
%     33           20            0.42876    0.42011      28.292        9.65        '1   2   4   6   8   9  10  13  14  15  19  20  21  23  24'             8    12     '{@covMaterniso, 3}'    'allPoints'                  1         '15*dim'           'meanConst'
%     51           17            0.44208    0.43416      33.458      9.2941        '3   4   6   7   8  11  15  16  17  18  19  21  22  23  24'             5    12     '{@covMaterniso, 3}'    'clustering'             0.999         '20*dim'           'meanConst'
%     62           18            0.42152    0.39957      30.125      15.389        '1   2   3   4   5   7  10  11  12  13  14  15  17  19  21  22'         6    12     '{@covMaterniso, 5}'    'nearestToPopulation'        1         '10*dim'           'meanConst'

parameterSets_3a5l(1).covFcn = '{@covMaterniso, 5}';
parameterSets_3a5l(1).trainsetType = 'nearest';
parameterSets_3a5l(1).trainRange = 1;
parameterSets_3a5l(1).trainsetSizeMax = '20*dim';
parameterSets_3a5l(1).meanFcn = 'meanConst';
parameterSets_3a5l(1).trainAlgorithm = 'fmincon';
parameterSets_3a5l(1).hyp.lik = -4.605170;
parameterSets_3a5l(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_3a5l(2).covFcn = '{@covMaterniso, 3}';
parameterSets_3a5l(2).trainsetType = 'allPoints';
parameterSets_3a5l(2).trainRange = 1;
parameterSets_3a5l(2).trainsetSizeMax = '10*dim';
parameterSets_3a5l(2).meanFcn = 'meanConst';
parameterSets_3a5l(2).trainAlgorithm = 'fmincon';
parameterSets_3a5l(2).hyp.lik = -4.605170;
parameterSets_3a5l(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_3a5l(3).covFcn = '{@covMaterniso, 5}';
parameterSets_3a5l(3).trainsetType = 'allPoints';
parameterSets_3a5l(3).trainRange = 1;
parameterSets_3a5l(3).trainsetSizeMax = '10*dim';
parameterSets_3a5l(3).meanFcn = 'meanConst';
parameterSets_3a5l(3).trainAlgorithm = 'fmincon';
parameterSets_3a5l(3).hyp.lik = -4.605170;
parameterSets_3a5l(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_3a5l(4).covFcn = '{@covMaterniso, 3}';
parameterSets_3a5l(4).trainsetType = 'allPoints';
parameterSets_3a5l(4).trainRange = 1;
parameterSets_3a5l(4).trainsetSizeMax = '15*dim';
parameterSets_3a5l(4).meanFcn = 'meanConst';
parameterSets_3a5l(4).trainAlgorithm = 'fmincon';
parameterSets_3a5l(4).hyp.lik = -4.605170;
parameterSets_3a5l(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_3a5l(5).covFcn = '{@covMaterniso, 3}';
parameterSets_3a5l(5).trainsetType = 'clustering';
parameterSets_3a5l(5).trainRange = 0.999;
parameterSets_3a5l(5).trainsetSizeMax = '20*dim';
parameterSets_3a5l(5).meanFcn = 'meanConst';
parameterSets_3a5l(5).trainAlgorithm = 'fmincon';
parameterSets_3a5l(5).hyp.lik = -4.605170;
parameterSets_3a5l(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_3a5l(6).covFcn = '{@covMaterniso, 5}';
parameterSets_3a5l(6).trainsetType = 'nearestToPopulation';
parameterSets_3a5l(6).trainRange = 1;
parameterSets_3a5l(6).trainsetSizeMax = '10*dim';
parameterSets_3a5l(6).meanFcn = 'meanConst';
parameterSets_3a5l(6).trainAlgorithm = 'fmincon';
parameterSets_3a5l(6).hyp.lik = -4.605170;
parameterSets_3a5l(6).hyp.cov = [-0.693147; 0.693147];

modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', {0.25},...
    'maxGenerationShiftForModelChoice', {2},...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_3a5l },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
