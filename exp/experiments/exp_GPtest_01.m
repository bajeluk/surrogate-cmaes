% exp_GPtest_01 model testing experiment -- Matlab part

% EXP_ID
exp_id = 'exp_GPtest_01';

% FUN/DIM/INST settings
if ~exist('func', 'var')
  func = (1:24);  end
if ~exist('dims', 'var')
  dims = [2, 5, 10];  end
if ~exist('instances', 'var')
  instances = [1:5 41:50];  end

% Maximal number of function evaluation per dimension to consider
maxEvals = 100;

% path settings
scratch = getenv('SCRATCHDIR');

experimentDir   = fullfile('exp', 'experiments', exp_id);
datasetFilename = fullfile(experimentDir, 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);
% defSetFile = fullfile(scratch, 'model', 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);
defModelFolder = fullfile(experimentDir, 'defData', ['defModel_', num2str(maxEvals), 'FE']);

% default model options
defModelOptions.useShift = false;
defModelOptions.predictionType = 'sd2';
defModelOptions.trainAlgorithm = 'fmincon';
defModelOptions.covFcn = '{@covMaterniso, 5}';
defModelOptions.normalizeY = true;
defModelOptions.hyp.lik = log(0.01);
defModelOptions.hyp.cov = log([0.5; 2]);

% test model options
modelSet = defModelOptions;
% none binning
modelSet.prediction = {'avgord'};
modelSet.binning = 'none';
modelOptions = combineFieldValues(modelSet);
% the rest of settings
modelSet.binning = {'logcluster', 'unipoints'};
modelSet.nBins = {'mu', 'lambda', '2*lambda'};
modelRest = combineFieldValues(modelSet);
% add rest to model options
modelOptions = [modelOptions; modelRest];

%% create testing dataset
% modelTestSets('exp_doubleEC_21_log', func, dims, maxEvals);
% load dataset
% ds_load = load(fullfile(experimentDir, 'defSet'));
% ds = ds_load;

%% test chosen models
modelFolders = testModels('ordgp', modelOptions, defSetFile, func, dims, false, experimentDir);
% add default folder (GpModel)
modelFolders = [defModelFolder; modelFolders];

%% compare results
modelStatistics(modelFolders, func, dims)
