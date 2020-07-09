% Script splitting generation samples from DTS-CMA-ES runs stored in file 
% DTS_meta_005.mat to settings validation (for RF), training, and testing 
% datasets (1:6:1). Splitting is
% performed at random on the model id level (8 models, 9th ADD excluded due
% to missing data in 20D).

%% load data
datasetFolder = fullfile('exp', 'experiments', 'dataset');
origFile = fullfile(datasetFolder, 'DTS_meta_005.mat');
dat = load(origFile);

%% Split data

% set seed for similar results
rng(42)

% exclude ADD covariance function
useModelSet = [1:6,8:9];
ds_new = dat.ds(:, :, :, useModelSet);
modelSettings_new = dat.modelSettings(useModelSet);
[nFunc, nDim, nInst, nMSet] = size(ds_new);
dsNewSize = [nFunc, nDim, nInst, nMSet];

% prepare resulting dataset cellarrays
ds_validation = cell(nFunc, nDim, nInst);
ds_train = cell(nFunc, nDim, nInst, 7);
ds_test = cell(nFunc, nDim, nInst);

% prepare matrix of all possible combinations of functions, dimensions, and
% instances
combMat = combvec(1:nFunc, 1:nDim, 1:nInst);

% generate random permutations of modelSettings
msPerm = cell2mat(arrayfun(@(x) randperm(nMSet, nMSet), (1:(nFunc*nDim*nInst))', 'UniformOutput', false));

% settings validation dataset
validCoor = [combMat', msPerm(:, 1)];
ds_validation(coor2ind(combMat', [nFunc, nDim, nInst])) = ds_new(coor2ind(validCoor, dsNewSize));

% training dataset
trainCoor = cell2mat(arrayfun(@(x) [combMat', msPerm(:, x)], 1:7, 'UniformOutput', false)');
traindsCoor = cell2mat(arrayfun(@(x) [combMat', x*ones(size(msPerm,1), 1)], 1:7, 'UniformOutput', false)');
ds_train(coor2ind(traindsCoor, [nFunc, nDim, nInst, 7])) = ds_new(coor2ind(trainCoor, dsNewSize));

% testing dataset
testCoor = [combMat', msPerm(:, 8)];
ds_test(coor2ind(combMat'), [nFunc, nDim, nInst]) = ds_new(coor2ind(testCoor, dsNewSize));

%% save resulting datasets

dat_new = dat;
dat_new.modelSettings = modelSettings_new;
dat_new.opts.note = 'Created using split_DTS_meta_data.m';

% validation dataset
dat_validation = dat_new;
dat_validation.ds = ds_validation;
dat_validation.opts.datasetName = 'DTS_meta_005_validation';
dat_validation.opts.datasetFile = fullfile(datasetFolder, 'DTS_meta_005_validation.mat');
save(dat_validation.opts.datasetFile, '-struct', 'dat_validation')

% training dataset
dat_train = dat_new;
dat_train.ds = ds_train;
dat_train.opts.datasetName = 'DTS_meta_005_train';
dat_train.opts.datasetFile = fullfile(datasetFolder, 'DTS_meta_005_train.mat');
save(dat_train.opts.datasetFile, '-struct', 'dat_train')

% testing dataset
dat_test = dat_new;
dat_test.ds = ds_test;
dat_test.opts.datasetName = 'DTS_meta_005_test';
dat_test.opts.datasetFile = fullfile(datasetFolder, 'DTS_meta_005_test.mat');
save(dat_test.opts.datasetFile, '-struct', 'dat_test')