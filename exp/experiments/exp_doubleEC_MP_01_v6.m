exp_id = 'exp_doubleEC_MP_01_v6';
exp_description = 'DTS-CMA-ES with 1 settings of ModelPools: w/o ARD, and w/o meanLinear, in 2D, Surrogate CMA-ES, fixed DTS 0.05 (merged code) with sd2, DTIterations={1}, PreSampleSize=0.75, criterion={sd2}, 2pop';

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

% w/o ARD and w/o meanLinear (total # = 96)

% rdeRankingOpts.maxRank = {25, 35};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };

% == Best set of settings #1 (id = 13)== 
%    best according to the avg_covrd_RDE
%
%   id    maxRank    colStat    minTrainedPerc    colSetCover                           f_weight                           includeARD    include5dim    includeMeanLinear    nSettings    avg_nCovered    avg_RDE    avg_covrd_RDE    avg_rank    avg_covrd_rank    avg_n_covrd_funs    n_ARD    n_Linear
%   __    _______    _______    ______________    ___________    ______________________________________________________    __________    ___________    _________________    _________    ____________    _______    _____________    ________    ______________    ________________    _____    ________
%
%   13    25         2           0.7              2              @(nCover,modelErrors)(nCover.^(4))./sum(modelErrors,2)    false         true           false                5                18          0.42571    0.40835          39.538      11.562              13.8              0        0       
%
%   settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                            covrd_funs                             covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn  
%   _________    __________    _______    _________    ________    __________    ___________________________________________________________    _________    ____________________    _____________________    __________    _______________    ___________
% 
%    1           18            0.43028    0.40625          41      12.889        '2   5   6   7   8   9  10  11  12  13  14  15  16  18  22'    10     8     '{@covSEiso}'           'nearest'                    1         '5*dim'            'meanConst'
%   34           19            0.41766    0.40003      39.083      9.2105        '1   2   5   6   8   9  10  11  12  14  18  20  21'            11     8     '{@covMaterniso, 5}'    'allPoints'                  1         '5*dim'            'meanConst'
%   37           20            0.42295    0.41509      33.458        10.4        '1   2   4   8   9  10  11  13  14  19  20  22  23  24'         9    11     '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst'
%   69           16            0.44208    0.43526       45.75       10.25        '3   4   6   7  15  16  17  18  19  21  22  23  24'             4    12     '{@covMaterniso, 3}'    'clustering'             0.999         '20*dim'           'meanConst'
%   90           17            0.41561    0.38512      38.396      15.059        '1   2   3   5   7   8  11  13  15  17  18  20  21  23'         5    12     '{@covMaterniso, 5}'    'nearestToPopulation'    0.999         '5*dim'            'meanConst'

parameterSets_13woAL(1).covFcn = '{@covSEiso}';
parameterSets_13woAL(1).trainsetType = 'nearest';
parameterSets_13woAL(1).trainRange = 1;
parameterSets_13woAL(1).trainsetSizeMax = '5*dim';
parameterSets_13woAL(1).meanFcn = 'meanConst';
parameterSets_13woAL(1).trainAlgorithm = 'fmincon';
parameterSets_13woAL(1).hyp.lik = -4.605170;
parameterSets_13woAL(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_13woAL(2).covFcn = '{@covMaterniso, 5}';
parameterSets_13woAL(2).trainsetType = 'allPoints';
parameterSets_13woAL(2).trainRange = 1;
parameterSets_13woAL(2).trainsetSizeMax = '5*dim';
parameterSets_13woAL(2).meanFcn = 'meanConst';
parameterSets_13woAL(2).trainAlgorithm = 'fmincon';
parameterSets_13woAL(2).hyp.lik = -4.605170;
parameterSets_13woAL(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_13woAL(3).covFcn = '{@covMaterniso, 3}';
parameterSets_13woAL(3).trainsetType = 'allPoints';
parameterSets_13woAL(3).trainRange = 1;
parameterSets_13woAL(3).trainsetSizeMax = '10*dim';
parameterSets_13woAL(3).meanFcn = 'meanConst';
parameterSets_13woAL(3).trainAlgorithm = 'fmincon';
parameterSets_13woAL(3).hyp.lik = -4.605170;
parameterSets_13woAL(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_13woAL(4).covFcn = '{@covMaterniso, 3}';
parameterSets_13woAL(4).trainsetType = 'clustering';
parameterSets_13woAL(4).trainRange = 0.999;
parameterSets_13woAL(4).trainsetSizeMax = '20*dim';
parameterSets_13woAL(4).meanFcn = 'meanConst';
parameterSets_13woAL(4).trainAlgorithm = 'fmincon';
parameterSets_13woAL(4).hyp.lik = -4.605170;
parameterSets_13woAL(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_13woAL(5).covFcn = '{@covMaterniso, 5}';
parameterSets_13woAL(5).trainsetType = 'nearestToPopulation';
parameterSets_13woAL(5).trainRange = 0.999;
parameterSets_13woAL(5).trainsetSizeMax = '5*dim';
parameterSets_13woAL(5).meanFcn = 'meanConst';
parameterSets_13woAL(5).trainAlgorithm = 'fmincon';
parameterSets_13woAL(5).hyp.lik = -4.605170;
parameterSets_13woAL(5).hyp.cov = [-0.693147; 0.693147];


% rdeRankingOpts.maxRank = {15, 25};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };


% == Best set of settings #1 (id = 1)== 
%    best according to the max_covrd_RDE
% 
%   id    maxRank    colStat    minTrainedPerc    colSetCover                           f_weight                           includeARD    include5dim    includeMeanLinear    nSettings    avg_nCovered    avg_RDE    avg_covrd_RDE    max_covrd_RDE    avg_rank    avg_covrd_rank    avg_n_covrd_funs    n_ARD    n_Linear
%   __    _______    _______    ______________    ___________    ______________________________________________________    __________    ___________    _________________    _________    ____________    _______    _____________    _____________    ________    ______________    ________________    _____    ________
%
%    1    15         2           0.7              1              @(nCover,modelErrors)(nCover.^(3))./sum(modelErrors,2)    false         true           false                10             10.5          0.42872    0.41333          0.43338          40.156      7.5266               9.4              0        0       
%
%   settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                      covrd_funs                       covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn  
%   _________    __________    _______    _________    ________    __________    _______________________________________________    _________    ____________________    _____________________    __________    _______________    ___________
%
%    5           10             0.4211    0.40966      40.729         6.3        '2   5   9  10  11  12  16  18  19'                9    1       '{@covSEiso}'           'allPoints'                  1         '5*dim'            'meanConst'
%   30           10            0.44526    0.43338      43.562        10.3        '3   4   6   7   8  17  18  21  24'                4    6       '{@covSEiso}'           'clustering'             0.999         '10*dim'           'meanConst'
%   34           15            0.41766    0.40077      39.083         6.2        '1   2   5   8   9  10  11  12  14  18  20  21'    9    6       '{@covMaterniso, 5}'    'allPoints'                  1         '5*dim'            'meanConst'
%   37           12            0.42295    0.41344      33.458      4.1667        '1   4   8   9  11  13  19  20  22  23'            7    5       '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst'
%   39            8             0.4305    0.40833      39.458       8.625        '1   8   9  10  13  20  23'                        3    5       '{@covMaterniso, 3}'    'allPoints'                  1         '20*dim'           'meanConst'
%   45           12            0.42876    0.42214      39.083      6.6667        '1   2   4   6   9  13  19  20  21  23  24'        8    4       '{@covMaterniso, 3}'    'allPoints'                  1         '15*dim'           'meanConst'
%   69           10            0.44208    0.43002       45.75         5.7        '3   4   7  15  16  18  19  22  23  24'            1    9       '{@covMaterniso, 3}'    'clustering'             0.999         '20*dim'           'meanConst'
%   70            9            0.44335    0.42324      46.812      8.4444        '3   4  15  16  18  21  22  24'                    1    8       '{@covMaterniso, 5}'    'clustering'             0.999         '20*dim'           'meanConst'
%   73            8             0.4193    0.39633      36.229         8.5        '2  11  12  13  14  18  23'                        4    4       '{@covSEiso}'           'nearestToPopulation'        1         '5*dim'            'meanConst'
%   82           11            0.41626    0.39603      37.396      10.364        '1   2   3   8  11  15  17  18  20  21  23'        3    8       '{@covMaterniso, 5}'    'nearestToPopulation'        1         '5*dim'            'meanConst'

parameterSets_1woAL(1).covFcn = '{@covSEiso}';
parameterSets_1woAL(1).trainsetType = 'allPoints';
parameterSets_1woAL(1).trainRange = 1;
parameterSets_1woAL(1).trainsetSizeMax = '5*dim';
parameterSets_1woAL(1).meanFcn = 'meanConst';
parameterSets_1woAL(1).trainAlgorithm = 'fmincon';
parameterSets_1woAL(1).hyp.lik = -4.605170;
parameterSets_1woAL(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(2).covFcn = '{@covSEiso}';
parameterSets_1woAL(2).trainsetType = 'clustering';
parameterSets_1woAL(2).trainRange = 0.999;
parameterSets_1woAL(2).trainsetSizeMax = '10*dim';
parameterSets_1woAL(2).meanFcn = 'meanConst';
parameterSets_1woAL(2).trainAlgorithm = 'fmincon';
parameterSets_1woAL(2).hyp.lik = -4.605170;
parameterSets_1woAL(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(3).covFcn = '{@covMaterniso, 5}';
parameterSets_1woAL(3).trainsetType = 'allPoints';
parameterSets_1woAL(3).trainRange = 1;
parameterSets_1woAL(3).trainsetSizeMax = '5*dim';
parameterSets_1woAL(3).meanFcn = 'meanConst';
parameterSets_1woAL(3).trainAlgorithm = 'fmincon';
parameterSets_1woAL(3).hyp.lik = -4.605170;
parameterSets_1woAL(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(4).covFcn = '{@covMaterniso, 3}';
parameterSets_1woAL(4).trainsetType = 'allPoints';
parameterSets_1woAL(4).trainRange = 1;
parameterSets_1woAL(4).trainsetSizeMax = '10*dim';
parameterSets_1woAL(4).meanFcn = 'meanConst';
parameterSets_1woAL(4).trainAlgorithm = 'fmincon';
parameterSets_1woAL(4).hyp.lik = -4.605170;
parameterSets_1woAL(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(5).covFcn = '{@covMaterniso, 3}';
parameterSets_1woAL(5).trainsetType = 'allPoints';
parameterSets_1woAL(5).trainRange = 1;
parameterSets_1woAL(5).trainsetSizeMax = '20*dim';
parameterSets_1woAL(5).meanFcn = 'meanConst';
parameterSets_1woAL(5).trainAlgorithm = 'fmincon';
parameterSets_1woAL(5).hyp.lik = -4.605170;
parameterSets_1woAL(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(6).covFcn = '{@covMaterniso, 3}';
parameterSets_1woAL(6).trainsetType = 'allPoints';
parameterSets_1woAL(6).trainRange = 1;
parameterSets_1woAL(6).trainsetSizeMax = '15*dim';
parameterSets_1woAL(6).meanFcn = 'meanConst';
parameterSets_1woAL(6).trainAlgorithm = 'fmincon';
parameterSets_1woAL(6).hyp.lik = -4.605170;
parameterSets_1woAL(6).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(7).covFcn = '{@covMaterniso, 3}';
parameterSets_1woAL(7).trainsetType = 'clustering';
parameterSets_1woAL(7).trainRange = 0.999;
parameterSets_1woAL(7).trainsetSizeMax = '20*dim';
parameterSets_1woAL(7).meanFcn = 'meanConst';
parameterSets_1woAL(7).trainAlgorithm = 'fmincon';
parameterSets_1woAL(7).hyp.lik = -4.605170;
parameterSets_1woAL(7).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(8).covFcn = '{@covMaterniso, 5}';
parameterSets_1woAL(8).trainsetType = 'clustering';
parameterSets_1woAL(8).trainRange = 0.999;
parameterSets_1woAL(8).trainsetSizeMax = '20*dim';
parameterSets_1woAL(8).meanFcn = 'meanConst';
parameterSets_1woAL(8).trainAlgorithm = 'fmincon';
parameterSets_1woAL(8).hyp.lik = -4.605170;
parameterSets_1woAL(8).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(9).covFcn = '{@covSEiso}';
parameterSets_1woAL(9).trainsetType = 'nearestToPopulation';
parameterSets_1woAL(9).trainRange = 1;
parameterSets_1woAL(9).trainsetSizeMax = '5*dim';
parameterSets_1woAL(9).meanFcn = 'meanConst';
parameterSets_1woAL(9).trainAlgorithm = 'fmincon';
parameterSets_1woAL(9).hyp.lik = -4.605170;
parameterSets_1woAL(9).hyp.cov = [-0.693147; 0.693147];
parameterSets_1woAL(10).covFcn = '{@covMaterniso, 5}';
parameterSets_1woAL(10).trainsetType = 'nearestToPopulation';
parameterSets_1woAL(10).trainRange = 1;
parameterSets_1woAL(10).trainsetSizeMax = '5*dim';
parameterSets_1woAL(10).meanFcn = 'meanConst';
parameterSets_1woAL(10).trainAlgorithm = 'fmincon';
parameterSets_1woAL(10).hyp.lik = -4.605170;
parameterSets_1woAL(10).hyp.cov = [-0.693147; 0.693147];


% rdeRankingOpts.maxRank = {20, 25};
% rdeRankingOpts.colStat = {2};
% rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% rdeRankingOpts.colSetCover = {1, 2};
% rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
%     @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };


% == Best set of settings #1 (id = 2)== 
%    best according to the max_covrd_RDE
% 
%   id    maxRank    colStat    minTrainedPerc    colSetCover                           f_weight                           includeARD    include5dim    includeMeanLinear    nSettings    avg_nCovered    avg_RDE    avg_covrd_RDE    max_covrd_RDE    avg_rank    avg_covrd_rank    avg_n_covrd_funs    n_ARD    n_Linear
%   __    _______    _______    ______________    ___________    ______________________________________________________    __________    ___________    _________________    _________    ____________    _______    _____________    _____________    ________    ______________    ________________    _____    ________
%
%    2    25         2           0.7              1              @(nCover,modelErrors)(nCover.^(3))./sum(modelErrors,2)    false         true           false                 6           16.167          0.42816    0.41259          0.43526          42.635      11.459            12.333              0        0       
%
%   settingNo    No_covered    avg_RDE    covrd_RDE    avg_rank    covrd_rank                          covrd_funs                           covrd_snp           covFcn               trainsetType         trainRange    trainsetSizeMax      meanFcn  
%   _________    __________    _______    _________    ________    __________    _______________________________________________________    _________    ____________________    _____________________    __________    _______________    ___________
%
%    5           18             0.4211    0.40909      40.729      12.833        '1   2   5   8   9  10  11  12  15  16  18  19  20  22'    11     7     '{@covSEiso}'           'allPoints'                  1         '5*dim'            'meanConst'
%   27            7            0.44959    0.43094      58.396          11        '5   6  10  12  14  16'                                     6     1     '{@covMaterniso, 5}'    'nearest'                    1         '20*dim'           'meanConst'
%   34           19            0.41766    0.40003      39.083      9.2105        '1   2   5   6   8   9  10  11  12  14  18  20  21'        11     8     '{@covMaterniso, 5}'    'allPoints'                  1         '5*dim'            'meanConst'
%   37           20            0.42295    0.41509      33.458        10.4        '1   2   4   8   9  10  11  13  14  19  20  22  23  24'     9    11     '{@covMaterniso, 3}'    'allPoints'                  1         '10*dim'           'meanConst'
%   69           16            0.44208    0.43526       45.75       10.25        '3   4   6   7  15  16  17  18  19  21  22  23  24'         4    12     '{@covMaterniso, 3}'    'clustering'             0.999         '20*dim'           'meanConst'
%   90           17            0.41561    0.38512      38.396      15.059        '1   2   3   5   7   8  11  13  15  17  18  20  21  23'     5    12     '{@covMaterniso, 5}'    'nearestToPopulation'    0.999         '5*dim'            'meanConst'

parameterSets_2woAL(1).covFcn = '{@covSEiso}';
parameterSets_2woAL(1).trainsetType = 'allPoints';
parameterSets_2woAL(1).trainRange = 1;
parameterSets_2woAL(1).trainsetSizeMax = '5*dim';
parameterSets_2woAL(1).meanFcn = 'meanConst';
parameterSets_2woAL(1).trainAlgorithm = 'fmincon';
parameterSets_2woAL(1).hyp.lik = -4.605170;
parameterSets_2woAL(1).hyp.cov = [-0.693147; 0.693147];
parameterSets_2woAL(2).covFcn = '{@covMaterniso, 5}';
parameterSets_2woAL(2).trainsetType = 'nearest';
parameterSets_2woAL(2).trainRange = 1;
parameterSets_2woAL(2).trainsetSizeMax = '20*dim';
parameterSets_2woAL(2).meanFcn = 'meanConst';
parameterSets_2woAL(2).trainAlgorithm = 'fmincon';
parameterSets_2woAL(2).hyp.lik = -4.605170;
parameterSets_2woAL(2).hyp.cov = [-0.693147; 0.693147];
parameterSets_2woAL(3).covFcn = '{@covMaterniso, 5}';
parameterSets_2woAL(3).trainsetType = 'allPoints';
parameterSets_2woAL(3).trainRange = 1;
parameterSets_2woAL(3).trainsetSizeMax = '5*dim';
parameterSets_2woAL(3).meanFcn = 'meanConst';
parameterSets_2woAL(3).trainAlgorithm = 'fmincon';
parameterSets_2woAL(3).hyp.lik = -4.605170;
parameterSets_2woAL(3).hyp.cov = [-0.693147; 0.693147];
parameterSets_2woAL(4).covFcn = '{@covMaterniso, 3}';
parameterSets_2woAL(4).trainsetType = 'allPoints';
parameterSets_2woAL(4).trainRange = 1;
parameterSets_2woAL(4).trainsetSizeMax = '10*dim';
parameterSets_2woAL(4).meanFcn = 'meanConst';
parameterSets_2woAL(4).trainAlgorithm = 'fmincon';
parameterSets_2woAL(4).hyp.lik = -4.605170;
parameterSets_2woAL(4).hyp.cov = [-0.693147; 0.693147];
parameterSets_2woAL(5).covFcn = '{@covMaterniso, 3}';
parameterSets_2woAL(5).trainsetType = 'clustering';
parameterSets_2woAL(5).trainRange = 0.999;
parameterSets_2woAL(5).trainsetSizeMax = '20*dim';
parameterSets_2woAL(5).meanFcn = 'meanConst';
parameterSets_2woAL(5).trainAlgorithm = 'fmincon';
parameterSets_2woAL(5).hyp.lik = -4.605170;
parameterSets_2woAL(5).hyp.cov = [-0.693147; 0.693147];
parameterSets_2woAL(6).covFcn = '{@covMaterniso, 5}';
parameterSets_2woAL(6).trainsetType = 'nearestToPopulation';
parameterSets_2woAL(6).trainRange = 0.999;
parameterSets_2woAL(6).trainsetSizeMax = '5*dim';
parameterSets_2woAL(6).meanFcn = 'meanConst';
parameterSets_2woAL(6).trainAlgorithm = 'fmincon';
parameterSets_2woAL(6).hyp.lik = -4.605170;
parameterSets_2woAL(6).hyp.cov = [-0.693147; 0.693147];

modelParams = { ...
    'retrainPeriod',      { 1 }, ...
    'bestModelSelection', { 'rdeAll' }, ...
    'historyLength',      { 7 }, ...
    'minTrainedModelsPercentilForModelChoice', { 0.25 },...
    'maxGenerationShiftForModelChoice', { 2 },...
    'predictionType',     { 'sd2' }, ...
    'useShift',           { false }, ...
    'normalizeY',         { true }, ...
    'parameterSets', { parameterSets_13woAL, parameterSets_1woAL, parameterSets_2woAL },...
    };


% CMA-ES parameters

cmaesParams = { ...
  'PopSize',            { '(8+floor(6*log(N)))' }, ...        %, '(8 + floor(6*log(N)))'};
  'Restarts',           { 50 }, ...
  'DispModulo',         { 0 }, ...
};

logDir = '/storage/plzen1/home/juranja3/public';
