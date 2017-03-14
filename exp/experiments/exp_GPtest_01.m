% exp_GPtest_01 model testing experiment -- Matlab part

% EXP_ID
exp_id  = 'exp_GPtest_01';
dataset = 'DTS_005'

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
datasetFilename = fullfile(experimentDir, 'dataset', [dataset, '.mat']);
% defSetFile = fullfile(scratch, 'model', 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);

% default model options
defModelOptions.useShift = false;
defModelOptions.predictionType = 'sd2';
defModelOptions.trainAlgorithm = 'fmincon';
defModelOptions.covFcn = '{@covMaterniso, 5}';
defModelOptions.normalizeY = true;
defModelOptions.hyp.lik = log(0.01);
defModelOptions.hyp.cov = log([0.5; 2]);
defModelOptions.covBounds = [ [-2;-2], [25;25] ];
defModelOptions.likBounds = log([1e-6, 10]);

models1 = defModelOptions;      % deep copy (!)
models1.trainsetType = { 'clustering', 'allPoints', 'nearest', 'nearestToPopulation' };
models1.trainRange   = { 0.99, 0.999 };
models1.trainsetSizeMax = { '5*dim', '10*dim', '15*dim', '20*dim' };
models1.meanFcn      = { 'meanConst', 'meanLinear' };
models1.covFcn       = { '{@covSEiso}', '{@covSEard}', ...
                         '{@covMaterniso, 5}', '{@covMaterniso, 3}' };

defModel_options = combineFieldValues(defModelOptions);
models1_options = combineFieldValues(models1);
modelOptions = [defModel_options; models1_options];

%% create testing dataset
% modelTestSets('exp_doubleEC_21_log', func, dims, maxEvals);
% load dataset
% ds_load = load(fullfile(experimentDir, 'defSet'));
% ds = ds_load;

%% test chosen models
modelFolders = testModels('ordgp', modelOptions, datasetFilename, func, dims, false, experimentDir);

%% compare results
modelStatistics(modelFolders, func, dims)
