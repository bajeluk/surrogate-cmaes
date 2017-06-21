exp_id = 'exp_BBComp_04';
exp_description = 'Surrogate CMA-ES settings for BBComp 1-OBJ-expensive task, DTS with preselection(1pt), DTIterations={1}, Valid.Gen.Period/PopSz=4/2, PreSampleSize=0.75, PoI criterion, 2pop';

% Surrogate manager parameters

surrogateParams.evoControl = 'doubletrained';
surrogateParams.observers = {'DTScreenStatistics', 'DTFileStatistics', 'NoneScreenStatistics', 'NoneFileStatistics', 'ECSaver'};
% surrogateParams.observers = {'NoneScreenStatistics'};
surrogateParams.modelType = 'gp';
surrogateParams.evoControlRestrictedParam = 0.05;
surrogateParams.evoControlTrainRange = '5+2*log2(dim)'; % '5*sqrt(dim)';             % will be multip. by sigma
surrogateParams.evoControlTrainNArchivePoints = '10*dim'; % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used

surrogateParams.updaterType = 'rankDiff';
surrogateParams.DTAdaptive_updateRate = 0.3;
surrogateParams.DTAdaptive_updateRateDown = 0.6;
surrogateParams.DTAdaptive_maxRatio = 0.05;
surrogateParams.DTAdaptive_minRatio = 0.05;
surrogateParams.DTAdaptive_defaultErr = 0.05;

surrogateParams.evoControlSwitchMode = 'none';
surrogateParams.evoControlSwitchTime = 3*24*3600; % 7*24*3600;
surrogateParams.evoControlMaxDoubleTrainIterations = 1;
surrogateParams.evoControlPreSampleSize = 0.75;
% surrogateParams.evoControlNBestPoints = [0.2 1.0];
% surrogateParams.evoControlValidationGenerationPeriod = 4;
surrogateParams.evoControlValidationPopSize = 0;
surrogateParams.evoControlOrigPointsRoundFcn = 'ceil'; % 'ceil', 'getProbNumber'

surrogateParams.evoControlIndividualExtension = [];    % will be multip. by lambda
surrogateParams.evoControlBestFromExtension = [];      % ratio of expanded popul.
surrogateParams.evoControlSampleRange = 1;             % will be multip. by sigma
surrogateParams.evoControlValidatePoints = [];

% Model parameters

surrogateParams.modelOpts.useShift        = false;
surrogateParams.modelOpts.predictionType  = 'poi';
surrogateParams.modelOpts.trainAlgorithm  = 'fmincon';
surrogateParams.modelOpts.covFcn          = '{@covMaterniso, 5}';
surrogateParams.modelOpts.normalizeY      = true;
surrogateParams.modelOpts.hyp.lik         = log(0.01);
surrogateParams.modelOpts.hyp.cov         = log([0.5; 2]);
surrogateParams.modelOpts.covBounds       = [ [-2;-2], [25;25] ];
surrogateParams.modelOpts.likBounds       = log([1e-6, 10]);

% EC Saver parameters
surrogateParams.maxArchSaveLen = 1e6;

% CMA-ES parameters

cmaesParams.PopSize = '(8 + floor(6*log(N)))';
cmaesParams.Restarts = 4;
cmaesParams.DispModulo = 50;
cmaesParams.EvalInitialX = true;
cmaesParams.EvalFinalMeanBeforeRestart = false;
cmaesParams.SaveVariables = 'on';
cmaesParams.SaveFilename = '[datapath filesep surrogateParams.exp_id ''_cmaesvars_'' surrogateParams.expFileID]';

% BOBYQA parameters

bobyqaParams.rho_beg = 0.3;
bobyqaParams.rho_end = 1e-7;
bobyqaParams.LBounds = 0;
bobyqaParams.UBounds = 1;

% BBCOMP client parameters

bbcompParams.libpath = 'exp/vendor/bbcomp/library/';
bbcompParams.libname = 'libbbcomp';
bbcompParams.libhfile = 'exp/vendor/bbcomp/client_matlab/bbcomplib.h';
bbcompParams.username = 'FILLME';
bbcompParams.password = 'FILLME';
bbcompParams.loghistory = 1; % set to 0 to disable history
bbcompParams.logfilepath = '[datapath filesep ''proxy_logs'']';
bbcompParams.trackname = 'BBComp2017-1OBJ-expensive';
bbcompParams.proxyHostname = 'localhost';
bbcompParams.proxyPort = 21000;
bbcompParams.proxyTimeout = 30;
bbcompParams.proxyConnectTimeout = 30;
bbcompParams.maxTrials = 1e2; % maximum no. of trials for network operations
bbcompParams.loginDelay = 2; % delay after each login retrial in seconds
bbcompParams.tryRecovery = false; % set to false to disable recovery
