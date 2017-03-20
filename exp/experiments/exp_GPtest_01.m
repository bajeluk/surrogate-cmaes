% exp_GPtest_01 model testing experiment -- Matlab part

% model training experiment settings
%
% EXPID (filename w/o extension of this script) | string
opts.exp_id     = 'exp_GPtest_01';
% dataset name (filename w/o extension of the dataset) | string
% or struct with the field '.ds' with 3D cell array with the data for {d, f, i}'s
opts.dataset    = 'DTS_005';
% type of model to test according to the ModelFactory | string
opts.modelType = 'gp';
% EXPPATH_SHORT
opts.exppath_short = fullfile('exp', 'experiments');
% statistics to compute
opts.statistics = { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde' };

% FUN/DIM/INST settings
if (~exist('func', 'var') || isempty(func))
  func = (1:24);  end
if (~exist('dims', 'var') || isempty(dims))
  dims = [2, 5, 10];  end
if (~exist('instances', 'var') || isempty(instances))
  instances = [1:5 41:50];  end

% Maximal number of function evaluation per dimension to consider
opts.maxEvals = 250;

% path settings
opts.scratch = getenv('SCRATCHDIR');

% other settings
opts.rewrite_results = false;

% directory with the dataset and results
opts.exppath = fullfile('exp', 'experiments', opts.exp_id);
% specifying the dataset -- expand the filename if dataset is string
if (ischar(opts.dataset))
  opts.dataset = fullfile(opts.exppath, 'dataset', [opts.dataset, '.mat']);
end
% defSetFile = fullfile(scratch, 'model', 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);

% default model options
defModelOptions.useShift        = false;
defModelOptions.predictionType  = 'sd2';
defModelOptions.trainAlgorithm  = 'fmincon';
defModelOptions.covFcn          = '{@covMaterniso, 5}';
defModelOptions.normalizeY      = true;
defModelOptions.hyp.lik         = log(0.01);
defModelOptions.hyp.cov         = log([0.5; 2]);
defModelOptions.covBounds       = [ [-2;-2], [25;25] ];
defModelOptions.likBounds       = log([1e-6, 10]);

% Full factorial design of the following parameters
models1 = defModelOptions;      % deep copy (!)
models1.trainsetType    = { 'nearestToPopulation', 'nearest', 'clustering', 'allPoints' };
models1.trainRange      = { 1.0, 0.999 };
models1.trainsetSizeMax = { '5*dim', '10*dim', '15*dim', '20*dim' };
models1.meanFcn         = { 'meanConst', 'meanLinear' };
models1.covFcn          = { '{@covSEiso}', '{@covSEard}', ...
                            '{@covMaterniso, 5}', '{@covMaterniso, 3}' };

defModel_options = combineFieldValues(defModelOptions);
models1_options  = combineFieldValues(models1);

% Combine default options and full factorial design
%modelOptions     = [defModel_options; models1_options];
modelOptions = models1_options;

%% create testing dataset
%ds = modelTestSets('exp_doubleEC_21_log15', func, dims, instances, opts);

fprintf('== Summary of the testing assignment ==\n');
fprintf('   # of models:  %d\n', length(modelOptions));
fprintf('   functions:    %s\n', num2str(func));
fprintf('   dimensions:   %s\n', num2str(dims));
fprintf('   instances:    %s\n', num2str(instances));
fprintf('=======================================\n');

%% test chosen models
modelFolders = testModels(modelOptions, opts, func, dims, instances);

%% load and calculate results
[rdeTable, mseTable, RDEs, MSEs] = modelStatistics(modelFolders, func, dims, instances);
