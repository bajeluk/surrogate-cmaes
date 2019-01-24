%% GECCO 2019 article plots
% Script for analysing GP covariance function dependence on data sampled 
% from DTS-CMA-ES run over the noiseless part of the BBOB framework.
% Script also creates graphs and tables.
% 
% Created for GECCO 2019 main track article.

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'gecco2019_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else

%%

% create tables from calculated features and save them

expfolder = fullfile('exp', 'experiments');
exp_id = 'exp_DTSmodels_meta_02';
exp_meta_output = fullfile(expfolder, exp_id, 'meta_output', 'exp_DTSmodels_meta_02_res.mat');
printScriptMess = false;

if ~isfile(exp_meta_output)
  if printScriptMess
    fprintf('Loading metalearing results\n')
  end
  [results, settings, resParams] = metaLearn_loadResults(fullfile(expfolder, exp_id), ...
    'ShowOutput', true, ...
    'SaveResults', exp_meta_output);
end

%%

% Load calculated data

if printScriptMess
  fprintf('Loading results\n')
end
if isfile(exp_meta_output)
  res = load(exp_meta_output, 'results', 'settings', 'resParams');
else
  error('%s does not exist', exp_meta_output)
end

results = res.results;
mfts = res.resParams.mfts;

nModel = numel(unique(results.model));
% create model labels
modelLabels = {'NN', 'SE', 'LIN', 'QUAD', ...
               'Matern', 'RQ', 'SE + QUAD', 'Gibbs'};             

%% RDE ranks 

if printScriptMess
  fprintf('Calculating RDE ranks\n')
end

rdeRankFile = fullfile('/tmp', 'normSumRDE.mat');
if ~exist('normSumRank', 'var') || ~isfile(rdeRankFile)
  % find unique data results
  [uniData, ia, ic] = unique(results(:, 1:5), 'rows');

  sumRank = zeros(nModel); % model x position
  for i = 1:size(ia, 1)
    if printScriptMess
      fprintf('%d/%d\n', i, size(ia, 1))
    end
    modelNum = results.model(ic == i);
    % use only results where 2 or more models are available
    if numel(modelNum) > 1
      values = results.rde(ic == i);
      % precise rank (from createRankingTable)
      [~, ~, tr] = unique(values);
      tr = tr';
      pr = arrayfun(@(x) sum(tr < x), tr) + 1;
      % add positions
      for m = 1:numel(modelNum)
        sumRank(modelNum(m), pr(m)) = sumRank(modelNum(m), pr(m)) + 1;
      end
    end
  end

  % normalize due to different numbers of comparisons (sometimes only two
  % models on one dataset are available)
  normSumRank = sumRank.*repmat(1./sum(sumRank, 2), 1, nModel);
  normSumRank = array2table(normSumRank, 'RowNames', modelLabels, ...
                'VariableNames', {'sum_1st', 'sum_2nd', 'sum_3rd', 'sum_4th', ...
                                  'sum_5th', 'sum_6th', 'sum_7th', 'sum_8th'});
  % save RDE rank result
  save(rdeRankFile, 'normSumRank')
elseif ~exist('normSumRank', 'var') && isfile(rdeRankFile)
  nsRank = load(rdeRankFile);
  normSumRank = nsRank.normSumRank;
end
disp(normSumRank)

%% Removing metafeatures
% Remove user-defined, constant, NaN, and linearly dependent metafeatures.

% Process metafeatures

if printScriptMess
  fprintf('Processing metafeatures\n')
end
% remove user-defined metafeatures

% list of excluded metafeatures
%  3: basic.lower_min
%  4: basic.lower_max
%  5: basic.upper_min
%  6: basic.upper_max
%  7: basic.objective_min
%  8: basic.objective_max
%  9: basic.blocks_min
% 10: basic.blocks_max
% 11: basic.cells_total
% 12: basic.cells_filled
% 13: basic.minimize_fun

mfts_sets = {'archive', 'train', 'traintest'};
rem_mfts = {};
for ms = 1:numel(mfts_sets)
  mss = [mfts_sets{ms}, '_'];
  rem_mfts = [rem_mfts, {...
      [mss, 'basic_lower_min'], ...
      [mss, 'basic_lower_max'], ...
      [mss, 'basic_upper_min'], ...
      [mss, 'basic_upper_max'], ...
      [mss, 'basic_objective_min'], ...
      [mss, 'basic_objective_max'], ...
      [mss, 'basic_blocks_min'], ...
      [mss, 'basic_blocks_max'], ...
      [mss, 'basic_cells_total'], ...
      [mss, 'basic_cells_filled'], ...
      [mss, 'basic_minimize_fun'], ...
      [mss, 'ela_metamodel_lin_simple_intercept'] ...
    }];
end
% also remove dimension of training set, because it cannot be determined by
% next steps
rem_mfts = [rem_mfts, 'train_basic_dim'];

fprintf('Removing user-defined metafeatures:\n')
for m = 1:numel(rem_mfts)
  mId = strcmp(mfts.Properties.VariableNames, rem_mfts{m});
  fprintf('%s\n', rem_mfts{m})
  mfts = mfts(:, ~mId);
end

%%

% Remove constant or NaN features
varMfts = varfun(@nanvar, mfts(:, 6:end), 'OutputFormat','uniform');
constId = varMfts == 0;
nanId = isnan(varMfts);

% constant
% add fun, dim, inst, id, and gen columns
constId = [false(1, 5), constId]; 
if any(constId)
  fprintf('Removing constant metafeatures:\n')
  constNumId = find(constId);
  for m = constNumId
    fprintf('%s\n', mfts.Properties.VariableNames{m})
  end
end

% NaN
% add fun, dim, inst, id, and gen columns
nanId = [false(1, 5), nanId]; 
if any(nanId)
  fprintf('Removing NaN metafeatures:\n')
  nanNumId = find(nanId);
  for m = nanNumId
    fprintf('%s\n', mfts.Properties.VariableNames{m})
  end
end
% actual removement
mfts(:, constId | nanId) = [];

%%

% Remove linearly dependent metafeatures

% tolerance for linear dependency
tol = eps; % eps;

mfts_only = mfts(:, 6:end);
mfts_arr = table2array(mfts_only);
nanInfMftsId = isnan(mfts_arr) + 2*isinf(mfts_arr);
% find columns where nan and Inf positions are identical
[~, nanInfDependent] = ismember(nanInfMftsId', nanInfMftsId', 'rows');
% cycle through identical positions to find linear dependencies of the
% remaining data
nanInfGroup = unique(nanInfDependent);
mfId = [];
for g = nanInfGroup'
  gid = g == nanInfDependent;
  % array consisting non-nan or inf rows of specific group
  subArray = mfts_arr(~nanInfMftsId(:, find(gid, 1, 'first')), gid);
  % exclude equal metafeatures
  [~, gEqId] = intersect(subArray', subArray', 'rows');
  % normalize columns
%   subArray = subArray./repmat(sum(subArray), size(subArray, 1), 1);
  % find linear independent columns of reduced group table
  [~, gDepId] = licols(subArray(:, gEqId), tol);
  gidNum = find(gid);
  % add independent metafeatures including NaNs and Infs
  mfId = [mfId; gidNum(gEqId(gDepId))];
end
% resulting table
mfts_indep = [mfts(:, 1:5), mfts_only(:, sort(mfId))];

nFeat = numel(mfId);
% print dependent metafeatures
if nFeat < numel(mfts_only.Properties.VariableNames)
  mfts_dep = mfts_only.Properties.VariableNames(...
    ~ismember(mfts_only.Properties.VariableNames, ...
              mfts_indep.Properties.VariableNames));
  fprintf('Removing dependent metafeatures:\n')
  for m = 1:numel(mfts_dep)
    fprintf('%s\n', mfts_dep{m})
  end
end

%%

% Cat model results with metafeatures

if printScriptMess
  fprintf('Concatenating model results and metafeatures\n')
end
err_name = 'rde';
errPenalty = 1;
% form table containing rde results from all models
for m = 1:nModel
  model_err = results(results.model == m, {'fun', 'dim', 'inst', 'gen', 'id', err_name});
  model_err.Properties.VariableNames{err_name} = sprintf('model%d_%s', m, err_name);
  if m == 1
    err_table = model_err;
  else
    err_table = outerjoin(err_table, model_err, 'MergeKeys', true);
  end
end
% replace missing errors and inf errors by penalty term
% cannot use fillmissing(err_table, 'constant', errPenalty) - not 
% implemented in R2015b
for m = 1:nModel
  col_name = sprintf('model%d_%s', m, err_name);
  err_table.(col_name)(isnan(err_table.(col_name))) = errPenalty;
end

% join metafeatures and error values of individual models
mfts_err = innerjoin(err_table, mfts_indep);

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Visual inspection
% Medians (thick lines) and quartiles (thin and dash-dotted lines) of GP 
% model RDE dependency on individual metafeatures.

if printScriptMess
  fprintf('Starting visual inspection\n')
end
nPointsToPlot = 200;
logBound = 5;
medianLineWidth = 1.8;
quartileLineWidth = 1;
medianLineStyle = '-';
quartileLineStyle = '-.';

kerColor = getAlgColors([1, 2, 3, 12, 10, 11, 8, 5]) / 255;

close all

% metafeaturePlot(mfts_err.train_dispersion_diff_mean_02, mfts_err(:, 6:13), ...
metafeaturePlot(mfts_err(:, 14:end), mfts_err(:, 6:13), ...
                'DataColor', kerColor, ...
                'DataNames', modelLabels, ...
                'MftsNames', mfts_err.Properties.VariableNames(14:end), ...
                'NValues', 200, ...
                'LogBound', 2, ...
                'QuantileRange', [0.05, 0.95], ...
                'MedianLW', 1.8, ...
                'QuartileLW', 1, ...
                'MedianLS', '-', ...
                'QuartileLS', '-.' ...
  );

%% Correlation analysis
% Spearman correlation coefficients of Gaussian process model prediction
% RDE and ELA features.

if printScriptMess
  fprintf('Starting correlation analysis\n')
end
err_corr = NaN(nFeat, nModel);
for m = 1:nModel
  model_err_name = sprintf('model%d_%s', m, err_name);
  for mf = 1:nFeat
    mfts_name = mfts_indep.Properties.VariableNames{mf+5};
    % remove metafeature NaN and Inf values
    nanOrInf = isnan(mfts_err.(mfts_name)) | isinf(mfts_err.(mfts_name));
    actual_mfts_err = mfts_err.(mfts_name)(~nanOrInf);
    actual_model_err = mfts_err.(model_err_name)(~nanOrInf);
    err_corr(mf, m) = corr(actual_mfts_err, actual_model_err, ...
                           'Type', 'Spearman');
  end
end

% create feature labels
featureLabels = mfts_indep.Properties.VariableNames(6:end);
featureLabels = cellfun(@(x) strrep(x, '_', '\_'), featureLabels, 'UniformOutput', false);

% draw coeffs
close all
% image settings
labelRot = 60;
sizeX = 20;
sizeY = 46;
% draw coeffs without colorbar
han{1} = figure('Units', 'centimeters', ...
                'Position', [1 1 sizeX sizeY], ...
                'PaperSize', [sizeX + 2, sizeY + 2]);
imagesc(err_corr);
colorbar
% axis square
ax = gca;
ax.YTick = 1:nFeat;
ax.YTickLabel = featureLabels;
ax.XTick = 1:nModel;
ax.XTickLabel = modelLabels;
ax.XTickLabelRotation = labelRot;

% print result to file
imageFolder = fullfile('test', 'local', 'images');
[~, ~] = mkdir(imageFolder);
resPdf = fullfile(imageFolder, 'gp_cov_meta_DTS_01_spearman.pdf');
print2pdf(han, resPdf, 1)

%%

% Decision tree analysis

if printScriptMess
  fprintf('Starting decision tree analysis\n')
end
% mfts array with penalty terms
mfts_arr_pen = real(table2array(mfts_err(:, 6+nModel:end)));

model_err = table2array(mfts_err(:, 5 + (1:nModel)));

% % case where we ignore that more models can have the same performance
% [~, best_model_id] = min(model_err, [], 2);
% 
% CT = fitctree(mfts_arr_pen, best_model_id);
% CT = CT.prune('Level', 80);

%% 

% Regression tree analysis

% if printScriptMess
%   fprintf('Starting regression tree analysis\n')
% end
% nFolds = 5;
% 
% RT = cell(1, nModel);
% cvErrRT = NaN(nModel, nFolds);
% cvErrRF_100 = NaN(nModel, nFolds);
% cvId = cvInd(size(mfts_arr_pen, 1), nFolds);
% for m = 1:nModel
%   for f = 1:nFolds
%     if printScriptMess
%       fprintf('Model %d Fold %d\n', m, f)
%     end
%     trainId = cvId ~= f;
%     testId = ~trainId;
%     fprintf('Tree model\n')
%     % regression tree
%     actRT = fitrtree(mfts_arr_pen(trainId, :), model_err(trainId, m));
%     errPred = actRT.predict(mfts_arr_pen(testId, :));
%     cvErrRT(m, f) = mseLossFunc(errPred, model_err(testId, m));
%     fprintf('Forest model\n')
%     % regression forest
%     tic
%     actRF = TreeBagger(100, mfts_arr_pen(trainId, :), model_err(trainId, m), ...
%                        'Method', 'regression');
%     errPred = actRF.predict(mfts_arr_pen(testId, :));
%     cvErrRF_100(m, f) = mseLossFunc(errPred, model_err(testId, m));
%     toc
%   end
%   if printScriptMess
%     fprintf('Model %d RT training\n', m)
%   end
% %   RT{m} = fitrtree(mfts_arr_pen, model_err(:, m));
% end

%% Multi-label classification

if printScriptMess
  fprintf('Starting multi-label classification\n')
end
% find which models have error smaller than 0.1 quantile 
good_model = model_err <= repmat(quantile(model_err, 0.00, 2), 1, 8);
% cases where all models are bad
good_model(model_err == errPenalty) = false;

uni_good_model = unique(good_model, 'rows');
% print number of different classes
fprintf('Results contain %d of different classes:\n', ...
        size(uni_good_model, 1))
for cl = 1:size(uni_good_model)
  fprintf('%d: %s\n', cl, ...
    printStructure(modelLabels(logical(uni_good_model(cl, :))), 'Format', 'values'))
end

% train classifier for each model separately
for m = 1:nModel
  if printScriptMess
    fprintf('Training classifier %d\n', m)
  end
%   tic
%   modelClassifier{m} = TreeBagger(100, mfts_arr_pen, good_model(:, m));
%   mlcPred(:, m) = modelClassifier{m}.predict(mfts_arr_pen);
%   modelDT{m} = fitctree(mfts_arr_pen, good_model(:, m));
%   mlcPredDT(:, m) = modelDT{m}.predict(mfts_arr_pen);
%   toc
end

%% One-way ANOVA

close all

if printScriptMess
  fprintf('Starting one-way ANOVA\n')
end

nAllErr = numel(model_err);
[anova_p, anova_tbl, anova_stats] = ...
  anova1(reshape(model_err, nAllErr, 1), ...
         reshape(repmat(modelLabels, size(model_err, 1), 1), nAllErr, 1));
[multcomp_c, multcomp_m, multcomp_h, multcomp_nms] = multcompare(anova_stats);

%%
% ANOVA table and multcompare figure visualising differences among model
% RDE on the whole dataset.

%% PCA

close all

pca_input = table2array(mfts_err(:, 6+nModel:end));
pca_input(any(isinf(pca_input) | isnan(pca_input), 2), :) = [];

fprintf('Number of data (rows) not containing NaN or Inf: %d/%d\n', ...
  size(pca_input, 1), size(mfts_err, 1))
% pca(pca_input)