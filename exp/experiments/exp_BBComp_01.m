exp_id = 'exp_BBComp_01';
exp_description = 'Surrogate CMA-ES settings for BBComp 1-OBJ task, DTS with preselection(1pt), DTIterations={2}, Valid.Gen.Period/PopSz=4/2, PreSampleSize=0.75, PoI criterion, 2pop';

% Surrogate manager parameters

surrogateParams.evoControl = 'doubletrained';
surrogateParams.observers = {'DTScreenStatistics', 'DTFileStatistics', 'ECSaver'};
% surrogateParams.observers = {'NoneScreenStatistics'};
surrogateParams.modelType = 'gp';
surrogateParams.evoControlRestrictedParam = 0.05;
surrogateParams.evoControlTrainRange = 10;             % will be multip. by sigma
surrogateParams.evoControlTrainNArchivePoints = '15*dim'; % will be myeval()'ed, 'nRequired', 'nEvaluated', 'lambda', 'dim' can be used

surrogateParams.updaterType = 'rankDiff';
surrogateParams.DTAdaptive_updateRate = 0.3;
surrogateParams.DTAdaptive_updateRateDown = 0.6;
surrogateParams.DTAdaptive_maxRatio = 0.9;
surrogateParams.DTAdaptive_minRatio = 0.04;
surrogateParams.DTAdaptive_lowErr = 0.15;
surrogateParams.DTAdaptive_highErr = 0.40;
surrogateParams.DTAdaptive_defaultErr = 0.10;

% surrogateParams.evoControlSwitchMode = 'none';
% surrogateParams.evoControlSwitchTime = 1; % 7*24*3600;
surrogateParams.evoControlMaxDoubleTrainIterations = 2;
surrogateParams.evoControlPreSampleSize = 0.75;
surrogateParams.evoControlNBestPoints = [0.2 1.0];
surrogateParams.evoControlValidationGenerationPeriod = 4;
surrogateParams.evoControlValidationPopSize = 2;
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
surrogateParams.maxArchSaveLen = 10^6;

% CMA-ES parameters

cmaesParams.PopSize = '(8 + floor(6*log(N)))';
cmaesParams.Restarts = 4;
cmaesParams.DispModulo = 50;
cmaesParams.EvalInitialX = false;
cmaesParams.EvalFinalMeanBeforeRestart = false;
cmaesParams.SaveVariables = 'on';
cmaesParams.SaveFilename = '[datapath filesep surrogateParams.exp_id ''_cmaesvars_'' surrogateParams.expFileID ''.mat'']';

% BBCOMP client

bbcompParams.libpath = 'exp/vendor/bbcomp/library/';
bbcompParams.libname = 'libbbcomp';
bbcompParams.libhfile = 'exp/vendor/bbcomp/client_matlab/bbcomplib.h';
bbcompParams.username = 'FILLME';
bbcompParams.password = 'FILLME';
bbcompParams.trackname = 'trial';
bbcompParams.proxyHostname = 'localhost';
bbcompParams.proxyPort = 20000;
bbcompParams.proxyTimeout = 10;
bbcompParams.proxyConnectTimeout = 10;
bbcompParams.maxTrials = 1e2; % maximum no. of trials for network operations
bbcompParams.loginDelay = 2; % delay after each login retrial in seconds
