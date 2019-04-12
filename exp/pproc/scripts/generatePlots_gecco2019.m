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

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2019paper');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);

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

clear actualFolder articleFolder

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
modelLabels = {'NN', 'SE', 'LIN', 'Q', ...
               'Mat', 'RQ', 'SE+Q', 'Gibbs'};

%% RDE ranks 

% if printScriptMess
%   fprintf('Calculating RDE ranks\n')
% end
% 
% rdeRankFile = fullfile('/tmp', 'normSumRDE.mat');
% if ~exist('normSumRank', 'var') || ~isfile(rdeRankFile)
%   % find unique data results
%   [uniData, ia, ic] = unique(results(:, 1:5), 'rows');
% 
%   sumRank = zeros(nModel); % model x position
%   for i = 1:size(ia, 1)
%     if printScriptMess
%       fprintf('%d/%d\n', i, size(ia, 1))
%     end
%     modelNum = results.model(ic == i);
%     % use only results where 2 or more models are available
%     if numel(modelNum) > 1
%       values = results.rde(ic == i);
%       % precise rank (from createRankingTable)
%       [~, ~, tr] = unique(values);
%       tr = tr';
%       pr = arrayfun(@(x) sum(tr < x), tr) + 1;
%       % add positions
%       for m = 1:numel(modelNum)
%         sumRank(modelNum(m), pr(m)) = sumRank(modelNum(m), pr(m)) + 1;
%       end
%     end
%   end
% 
%   % normalize due to different numbers of comparisons (sometimes only two
%   % models on one dataset are available)
%   normSumRank = sumRank.*repmat(1./sum(sumRank, 2), 1, nModel);
%   normSumRank = array2table(normSumRank, 'RowNames', modelLabels, ...
%                 'VariableNames', {'sum_1st', 'sum_2nd', 'sum_3rd', 'sum_4th', ...
%                                   'sum_5th', 'sum_6th', 'sum_7th', 'sum_8th'});
%   % save RDE rank result
%   save(rdeRankFile, 'normSumRank')
% elseif ~exist('normSumRank', 'var') && isfile(rdeRankFile)
%   nsRank = load(rdeRankFile);
%   normSumRank = nsRank.normSumRank;
% end
% disp(normSumRank)

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
% also remove only train and trainset metafeatures identical for all mfts
% sets
rem2_mfts = {};
for ms = 2:numel(mfts_sets)
  mss = [mfts_sets{ms}, '_'];
  rem2_mfts = [rem2_mfts, {...
      [mss, 'basic_dim'], ...
      [mss, 'cmaes_cma_evopath_c_norm'], ...
      [mss, 'cmaes_cma_evopath_s_norm'], ...
      [mss, 'cmaes_cma_generation'], ...
      [mss, 'cmaes_cma_restart'], ...
      [mss, 'cmaes_cma_step_size'] ...
    }];
end
rem_mfts = [rem_mfts, rem2_mfts];

fprintf('Removing user-defined metafeatures:\n')
for m = 1:numel(rem_mfts)
  mId = strcmp(mfts.Properties.VariableNames, rem_mfts{m});
  fprintf('%s\n', rem_mfts{m})
  mfts = mfts(:, ~mId);
end

% rename archive identical columns (see above)
mfts.Properties.VariableNames{'archive_basic_dim'} = 'dimension';
mfts.Properties.VariableNames{'archive_cmaes_cma_evopath_c_norm'} = ...
                                      'cmaes_evopath_c_norm';
mfts.Properties.VariableNames{'archive_cmaes_cma_evopath_s_norm'} = ...
                                      'cmaes_evopath_s_norm';
mfts.Properties.VariableNames{'archive_cmaes_cma_generation'} = ...
                                      'cmaes_generation';
mfts.Properties.VariableNames{'archive_cmaes_cma_restart'} = ...
                                      'cmaes_restart';
mfts.Properties.VariableNames{'archive_cmaes_cma_step_size'} = ...
                                      'cmaes_step_size';
% rename observations columns
mfts.Properties.VariableNames{'archive_basic_observations'} = 'archive_observations';
mfts.Properties.VariableNames{'train_basic_observations'}   = 'train_observations';
mfts.Properties.VariableNames{'traintest_basic_observations'} = 'traintest_observations';

% create mfts table, where only user defined metafeatures are not present
full_mfts = mfts;

clear rem_mfts rem2_mfts m mId mss

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

clear varMfts constId nanId constNumId nanNumId m

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

clear tol mfts_only mfts_arr nanInfMftsId nanInfDependent nanInfGroup mfId ...
      g gid subArray gEqId gDepId gidNum nFeat mfts_dep m

%%

% Cat model results with metafeatures

if printScriptMess
  fprintf('Concatenating model results and metafeatures\n')
end
err_name = {'rdeValid', 'rde'};
errPenalty = 1;
modelLabels_friend = {'NN', 'SE', 'LIN', 'QUAD', ...
                     'Matern', 'RQ', 'SE_QUAD', 'Gibbs'};
full_mfts_err = full_mfts;

for e = 1:numel(err_name)
  % form table containing rde results from all models
  for m = 1:nModel
    model_err = results(results.model == m, {'fun', 'dim', 'inst', 'gen', 'id', err_name{e}});
    model_err.Properties.VariableNames{err_name{e}} = ...
      sprintf('%s_%s', err_name{e}, modelLabels_friend{m});
    if m == 1
      err_table = model_err;
    else
      err_table = outerjoin(err_table, model_err, 'MergeKeys', true);
    end
  end
  % join all metafeatures and error values of individual models without
  % replacing missing errors
  full_mfts_err = innerjoin(err_table, full_mfts_err);
  % model error column names
%   modelLabels_err = cellfun(@(x) [err_name, '_', x], modelLabels_friend, ...
%                             'UniformOutput', false);
%   full_mfts_err.Properties.VariableNames = ...
%     [modelLabels_err, full_mfts_err.Properties.VariableNames(9:end)];
end
% create extra full_mfts_err for visualisation
full_mfts_vis = full_mfts_err;
% remove identification columns (fun, dim, inst, gen, model)
full_mfts_err = full_mfts_err(:, [1, 3, 5, 6:end]);

% replace missing errors and inf errors by penalty term
% cannot use fillmissing(err_table, 'constant', errPenalty) - not 
% implemented in R2015b
for m = 1:nModel
  col_name = sprintf('%s_%s', err_name{e}, modelLabels_friend{m});
  err_table.(col_name)(isnan(err_table.(col_name))) = errPenalty;
end

% join metafeatures and error values of individual models
mfts_err = innerjoin(err_table, mfts_indep);
full_mfts_vis = innerjoin(err_table, full_mfts_vis(:, [1:5, 22:end]));

clear err_name errPenalty modelLabels_friend full_mfts full_mfts_err e m ...
      model_err err_table col_name

%%
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

nExtF = 55;
mfts_order = [14, 16:18, 20:21, ... identical
              15, 80, 138, ... observations
              19, 78, 136, ... cma_mean_dist
              22, 79, 137, ... cma_lik
              reshape(22 + repmat([0,58,116]', 1, nExtF) ...
                + [(1:nExtF); (1:nExtF); (1:nExtF)], 1, 3*nExtF)];

close all
 
% metafeaturePlot(full_mfts_vis(:, mfts_order), full_mfts_vis(:, 6:13), ...
%                 'DataColor', kerColor, ...
%                 'DataNames', modelLabels, ...
%                 'MftsNames', full_mfts_vis.Properties.VariableNames(mfts_order), ...
%                 'NValues', 200, ...
%                 'LogBound', 2, ...
%                 'QuantileRange', [0.05, 0.95], ...
%                 'MedianLW', 1.8, ...
%                 'QuartileLW', 1, ...
%                 'MedianLS', '-', ...
%                 'QuartileLS', '-.' ...
%   );

% print archive_cma_lik
mftId = 22;
mftName = ' ';
pdfNames = {fullfile(plotResultsFolder, 'archive_cma_lik.pdf')};

han = metafeaturePlot(table2array(full_mfts_vis(:, 22)), full_mfts_vis(:, 6:13), ...
                'DataColor', kerColor, ...
                'DataNames', modelLabels, ...
                'MftsNames', mftName, ...
                'NValues', 200, ...
                'LogBound', 2, ...
                'QuantileRange', [0.05, 0.95], ...
                'MedianLW', 1.8, ...
                'QuartileLW', 1, ...
                'MedianLS', '-', ...
                'QuartileLS', '-.' ...
  );

print2pdf(han, pdfNames, 1)

%% Correlation analysis
% Spearman correlation coefficients of Gaussian process model prediction
% RDE and ELA features.

% if printScriptMess
%   fprintf('Starting correlation analysis\n')
% end
% err_corr = NaN(nFeat, nModel);
% for m = 1:nModel
%   model_err_name = sprintf('model%d_%s', m, err_name{e});
%   for mf = 1:nFeat
%     mfts_name = mfts_indep.Properties.VariableNames{mf+5};
%     % remove metafeature NaN and Inf values
%     nanOrInf = isnan(mfts_err.(mfts_name)) | isinf(mfts_err.(mfts_name));
%     actual_mfts_err = mfts_err.(mfts_name)(~nanOrInf);
%     actual_model_err = mfts_err.(model_err_name)(~nanOrInf);
%     err_corr(mf, m) = corr(actual_mfts_err, actual_model_err, ...
%                            'Type', 'Spearman');
%   end
% end
% 
% % create feature labels
% featureLabels = mfts_indep.Properties.VariableNames(6:end);
% featureLabels = cellfun(@(x) strrep(x, '_', '\_'), featureLabels, 'UniformOutput', false);
% 
% % draw coeffs
% close all
% % image settings
% labelRot = 60;
% sizeX = 20;
% sizeY = 46;
% % draw coeffs without colorbar
% han{1} = figure('Units', 'centimeters', ...
%                 'Position', [1 1 sizeX sizeY], ...
%                 'PaperSize', [sizeX + 2, sizeY + 2]);
% imagesc(err_corr);
% colorbar
% % axis square
% ax = gca;
% ax.YTick = 1:nFeat;
% ax.YTickLabel = featureLabels;
% ax.XTick = 1:nModel;
% ax.XTickLabel = modelLabels;
% ax.XTickLabelRotation = labelRot;
% 
% % print result to file
% imageFolder = fullfile('test', 'local', 'images');
% [~, ~] = mkdir(imageFolder);
% resPdf = fullfile(imageFolder, 'gp_cov_meta_DTS_01_spearman.pdf');
% print2pdf(han, resPdf, 1)

%%

% Decision tree analysis

% if printScriptMess
%   fprintf('Starting decision tree analysis\n')
% end
% % mfts array with penalty terms
% mfts_arr_pen = real(table2array(mfts_err(:, 6+nModel:end)));
% 
% model_err = table2array(mfts_err(:, 5 + (1:nModel)));

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

% if printScriptMess
%   fprintf('Starting multi-label classification\n')
% end
% % find which models have error smaller than 0.1 quantile 
% good_model = model_err <= repmat(quantile(model_err, 0.00, 2), 1, 8);
% % cases where all models are bad
% good_model(model_err == errPenalty) = false;
% 
% uni_good_model = unique(good_model, 'rows');
% % print number of different classes
% fprintf('Results contain %d of different classes:\n', ...
%         size(uni_good_model, 1))
% for cl = 1:size(uni_good_model)
%   fprintf('%d: %s\n', cl, ...
%     printStructure(modelLabels(logical(uni_good_model(cl, :))), 'Format', 'values'))
% end
% 
% % train classifier for each model separately
% for m = 1:nModel
%   if printScriptMess
%     fprintf('Training classifier %d\n', m)
%   end
%   tic
%   modelClassifier{m} = TreeBagger(100, mfts_arr_pen, good_model(:, m));
%   mlcPred(:, m) = modelClassifier{m}.predict(mfts_arr_pen);
%   modelDT{m} = fitctree(mfts_arr_pen, good_model(:, m));
%   mlcPredDT(:, m) = modelDT{m}.predict(mfts_arr_pen);
%   toc
% end

%% One-way ANOVA

% close all
% 
% if printScriptMess
%   fprintf('Starting one-way ANOVA\n')
% end
% 
% nAllErr = numel(model_err);
% [anova_p, anova_tbl, anova_stats] = ...
%   anova1(reshape(model_err, nAllErr, 1), ...
%          reshape(repmat(modelLabels, size(model_err, 1), 1), nAllErr, 1));
% [multcomp_c, multcomp_m, multcomp_h, multcomp_nms] = multcompare(anova_stats);

%%
% ANOVA table and multcompare figure visualising differences among model
% RDE on the whole dataset.

%% PCA

% close all
% 
% pca_input = table2array(mfts_err(:, 6+nModel:end));
% pca_input(any(isinf(pca_input) | isnan(pca_input), 2), :) = [];
% 
% fprintf('Number of data (rows) not containing NaN or Inf: %d/%d\n', ...
%   size(pca_input, 1), size(mfts_err, 1))
% pca(pca_input)

%% Kolmogorov-Smirnov test

% Martin's code:
% UsableTable = full_mfts_err(~all(isnan(table2array(full_mfts_err(:, 1:8))), 2), :);
% ColumnNames = full_mfts_err.Properties.VariableNames;
% MetafeatureNames = cell(4, 58);
% MetafeatureNames(1, 1:6) = ColumnNames([9, 11:13, 15, 16]);
% MetafeatureNames(2:4, :) = ColumnNames([14, 17, 10, 18:72; 73:130; 131:188]);
% ReorderedTable = UsableTable(:, [1:9, 11:13, 15, 16, 14, 17, 10, 18:188]);
% Best=zeros(size(ReorderedTable,1),8);
% for r=1:size(ReorderedTable,1)
%   for c=1:8
%     if ~isnan(ReorderedTable(r,c))&&ReorderedTable(r,c)==min(ReorderedTable(r,~isnan(ReorderedTable(r,1:8))))
%       Best(r,c)=true;
%     end
%   end
% end
% for c=1:6
%   Distributions{1,c}=ReorderedTable(~isnan(ReorderedTable(:,c+8)),c+8);
% end
% for r=2:4
%   for c=1:58
%     Distributions{r,c}=ReorderedTable(~isnan(ReorderedTable(:,c+14+(r-2)*58)),c+14+(r-2)*58);
%   end
% end
% UncorrectedSignificance(1,7:58,1:8)=Inf;
% UncorrectedSignificance(4,54:58,1:8)=Inf;
% for m=1:8
%   for c=1:6
%     [~, UncorrectedSignificance(1,c,m)] = ...
%       kstest2(ReorderedTable(~isnan(ReorderedTable(:, c+8)) & ...
%                              Best(:, m), c+8), ...
%               Distributions{1, c});
%   end
%   for r=2:4
%     for c=1:53+5*sign(4-r)
%       [~, UncorrectedSignificance(r,c,m)] = ...
%         kstest2(ReorderedTable(~isnan(ReorderedTable(:, c+14+(r-2)*58)) & ...
%                                Best(:, m), c+14+(r-2)*58), ...
%                 Distributions{r,c});
%     end
%   end
% end
% VectorizedUncorrectedSignificance = reshape(UncorrectedSignificance, 4*58*8, 1);
% VectorizedCorrectedSignificance(1:1856, 1) = Inf;
% VectorizedCorrectedSignificance(isfinite(VectorizedUncorrectedSignificance)) = ...
%   bonf_holm(VectorizedUncorrectedSignificance(isfinite(VectorizedUncorrectedSignificance)));
% CorrectedSignificance = reshape(VectorizedCorrectedSignificance, 4, 58, 8);
% Significant05 = CorrectedSignificance < 0.05;

% load Martin's results
ks_res_file = fullfile(expfolder, exp_id, 'kstest_meta.mat');
ks_res = load(ks_res_file, 'CorrectedSignificance', 'ReorderedTable', ...
                           'MetafeatureNames', 'Best', 'Distributions');
%%
% KS test table

% significance table
ks_sign = ks_res.CorrectedSignificance;
ks_sign(isinf(ks_sign)) = NaN;
% metafeatures identical accross origin
nIdentical = sum(~isnan(ks_sign(1, :, 1)));
ks_sign = [[ks_sign(1, 1:nIdentical, :); NaN(2, nIdentical, nModel)], ...
           ks_sign(2:4, :, :)];
ks_sign_2 = [];
for m = 1:nModel
  ks_sign_2 = [ks_sign_2, ks_sign(:, :, m)'];
end
% replace p-values greater than one by one
ks_sign_2(ks_sign_2 > 1) = 1;

% metafeature names

% ks table setting
ksTabSet.ResultFile = fullfile(tableFolder, 'ksTable.tex');
ksTabSet.ColGroups = {'NN', 'SE', 'LIN', 'Q', 'Mat', 'RQ', 'SE + Q', 'Gibbs'};
ksTabSet.ColGroupNum = 3*ones(1, nModel);
ksTabSet.ColValName = '$\sampleset$';
ksTabSet.ColNames = repmat(...
                    {'$\archive$', '$\trainset$', '$\trainpredset$'}, ...
                    1, nModel);
featNames = { ...
    'dimension', ...
    ... cmaes
    'generation', ...
    'step_size', ...
    'restart', ...
    'evopath_c_norm', ...
    'evopath_s_norm', ...
    ...
    'cma_mean_dist', ...
    'cma_lik', ...
    ... 
    'observations', ...
    ... dispersion
    'ratio_mean_02', ...
    'ratio_median_02', ...
    'diff_mean_02', ...
    'diff_median_02', ...
    'ratio_mean_05', ...
    'ratio_median_05', ...
    'diff_mean_05', ...
    'diff_median_05', ...
    'ratio_mean_10', ...
    'ratio_median_10', ...
    'diff_mean_10', ...
    'diff_median_10', ...
    'ratio_mean_25', ...
    'ratio_median_25', ...
    'diff_mean_25', ...
    'diff_median_25', ...
    ... ela_distribution
    'skewness', ...
    'kurtosis', ...
    'num_of_peaks', ...
    ... ela_levelset
    'mmce_lda_10', ...
    'mmce_qda_10', ...
    'mmce_mda_10', ...
    'lda_qda_10', ...
    'lda_mda_10', ...
    'qda_mda_10', ...
    'mmce_lda_25', ...
    'mmce_qda_25', ...
    'mmce_mda_25', ...
    'lda_qda_25', ...
    'lda_mda_25', ...
    'qda_mda_25', ...
    'mmce_lda_50', ...
    'mmce_qda_50', ...
    'mmce_mda_50', ...
    'lda_qda_50', ...
    'lda_mda_50', ...
    'qda_mda_50', ...
    ... ela_metamodel
    'l_s_adj_r2', ... 'lin_simple_adj_r2', ...
    'l_s_c_min', ... 'lin_simple_coef_min', ...
    'l_s_c_max', ... 'lin_simple_coef_max', ...
    'l_s_c_max_min', ... 'lin_simple_coef_max_by_min', ...
    'l_w_i_adj_r2', ... 'lin_w_interact_adj_r2', ...
    'q_s_adj_r2', ... 'quad_simple_adj_r2', ...
    'q_s_cond', ... 'quad_simple_cond', ...
    'q_w_i_adj_r2', ... 'quad_w_interact_adj_r2', ...
    ... infocontent
    'h_max', ...
    'eps_s', ...
    'eps_max', ...
    'm0', ...
    'eps_ratio', ...
    ... nearest_better
    'nb_std_ratio', ...
    'nb_mean_ratio', ...
    'nb_cor', ...
    'dist_ratio', ...
    'nb_fitness_cor', ...
  };
% change _ to \_
tableFeatNames = cellfun(@(x) ['\\texttt{', strrep(x, '_', '\\_'), '}'], ...
                            featNames, ...
                            'UniformOutput', false);
featNonId = [nIdentical + [3,1,2], nIdentical+4 : numel(tableFeatNames)];
ksTabSet.RowNames = tableFeatNames(featNonId);
ksTabSet.RowGroups = {'\npt', '$\featCMA$', '$\featDisp$', '$\featyDis$', ...
                      '$\featLevel$', '$\featMM$', '$\featInfo$', '$\featNBC$', ...
                      };
ksTabSet.RowGroupNum = [1, 2, 16, 3, 18, 8, 5, 5];
ksTabSet.RowValName = '$\featset{}$';
ksTabSet.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality ',...
        'probability distributions of individual features on all data ', ...
        'and on those data on which a particular covariance function scored best. ', ...
        'P-values denote KS test results ', ...
        'rejecting the equality of both distributions with the Holm ', ...
        'correction at the family-wise significance level $\alpha=0.05$, ', ...
        'otherwise, p-values are not shown. ', ...
        '''X'' values in $\trainpredset$ column denote features impossible to calculate with points without fitness values. ', ...
        'Zeros indicate p-values below the smallest double precision number. ', ...
        '$\featMM$ notation: \texttt{l} -- \texttt{lin}, \texttt{q} -- \texttt{quad}, \texttt{s} -- \texttt{simple}, ', ...
        '\texttt{i} -- \texttt{interact}, \texttt{c} -- \texttt{coef}.'];
                          
% print ks table
prtSignifTable(ks_sign_2(featNonId, :), ksTabSet)

% settings for identical ks table
ksTabSetInd.ColGroups = '';
ksTabSetInd.ColGroupNum = [];
ksTabSetInd.ColNames = {'NN', 'SE', 'LIN', 'Q', 'Mat', 'RQ', 'SE+Q', 'Gibbs'};
ksTabSetInd.RowNames = tableFeatNames(1:nIdentical);
ksTabSetInd.RowGroups = {'\dm', '$\featCMA$'};
ksTabSetInd.RowGroupNum = [1, 5];
ksTabSetInd.ResultFile = fullfile(tableFolder, 'ksIndTable.tex');
ksTabSetInd.TableWidth = '\columnwidth';
ksTabSetInd.OneColumn = true;
ksTabSetInd.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality ',...
        'probability distributions of individual features on all data ', ...
        'and on those data on which a particular covariance function scored best. ', ...
        'Features in table pertain to all types of sample sets. ', ...
        'P-values denote KS test results ', ...
        'rejecting the equality of both distributions with the Holm ', ...
        'correction at the family-wise significance level $\alpha=0.05$, ', ...
        'otherwise, p-values are not shown. ', ...
        'Zeros indicate p-values below the smallest double precision number. '];

% print ks table for feature origin identical metafeatures
prtSignifTable(ks_sign_2(1:nIdentical, 1:3:end), ksTabSetInd)

%% 
% KS test image
close all

% image settings
sizeY = 17;
sizeX = 34;
labelRot = 60;
plotNames = {fullfile(plotResultsFolder, 'ks_fig.pdf')};

% change _ to \_
imgFeatNames = cellfun(@(x) strrep(x, '_', '\_'), ...
                            featNames, ...
                            'UniformOutput', false);

% draw coeffs without colorbar
han{1} = figure('Units', 'centimeters', ...
                'Position', [1 1 sizeX sizeY], ...
                'PaperSize', [sizeX + 2, sizeY + 2]);
% draw image
hold on
imagesc(-log(ks_sign_2(featNonId, :)'))
colorbar

% axis square
ax = gca;
ax.XTick = 1:numel(featNonId);
ax.XTickLabel = imgFeatNames(featNonId);
ax.YTick = 1:(3*nModel);
ax.YTickLabel = [modelLabels, modelLabels, modelLabels];
ax.XTickLabelRotation = labelRot;

hold off
print2pdf(han, plotNames, 1)

%% Distribution plots

close all
q_bound = [0.05, 0.95];

for r = 2 % 2:4
  for c = 20% 1:53 + 5*sign(4-r)
    distr_all = ks_res.Distributions{r,c};
    q_d_all = quantile(distr_all, q_bound);
    
    d_all_show = distr_all(distr_all > q_d_all(1) & distr_all < q_d_all(2));
    % x values for plot
    x_val = linspace(min(d_all_show), max(d_all_show));
%     x_val = logspace(log10(-min(d_all_show)), log10(-max(d_all_show)));
    pdca_all = fitdist(d_all_show, 'Kernel');    
    all_pdf = pdf(pdca_all, x_val);
    
          han = figure('Units', 'centimeters', ...
                   'Position', [1, 1, 16, 20], ...
                   'PaperSize', [16, 20]);
    
    % model loop
    for m = 1:nModel
      % values of covariance and sample set for distribution
      distr_cov = ks_res.ReorderedTable(~isnan(ks_res.ReorderedTable(:, c+14+(r-2)*58)) & ...
                                   ks_res.Best(:, m), c+14+(r-2)*58);
    
      q_d_cov = quantile(distr_cov, q_bound);  
      % range
      d_cov_show = distr_cov(distr_cov > q_d_cov(1) & distr_cov < q_d_cov(2));
    
      % fit probability distribution
      pdca_cov = fitdist(d_cov_show, 'Kernel');    
      cov_pdf = pdf(pdca_cov, x_val);
    
%       han(m) = figure('PaperSize', [14, 12]);
      subplot(nModel/2, 2, m) 
      area(x_val, all_pdf, 'LineWidth', 2, ...
                           'FaceColor', 'r', ...
                           'EdgeColor', 'r', ...
                           'FaceAlpha', 0.2 ...
          )
      hold on
      area(x_val, cov_pdf, 'LineWidth', 2, ...
                           'FaceColor', 'b', ...
                           'EdgeColor', 'b', ...
                           'FaceAlpha', 0.2 ...
          )
      gca_act = gca;
      axis([min(d_all_show), max(d_all_show), 0, 0.25])
      if mod(m, 2) == 1
        ylabel('PDF')
      end
%       title([ks_res.MetafeatureNames{r, c}, ' for ', modelLabels{m}])
      title(modelLabels{m})
      legend({'all', modelLabels{m}})
      hold off
    end
%     h_cov = histfit(distr_cov(distr_cov > q_d_cov(1) & distr_cov < q_d_cov(2)), ...
%                     100, 'kernel');
%     h_all = histfit(distr_all(distr_all > q_d_all(1) & distr_all < q_d_all(2)), ...
%                     100, 'kernel');
%     figure(4)
    
  end
end

% distrFigNames = cellfun(@(x) fullfile(plotResultsFolder, ['skew_', x, '.pdf']), ...
%                         modelLabels, 'UniformOutput', false);
distrFigNames = {fullfile(plotResultsFolder, 'archive_skewness.pdf')};
print2pdf(han, distrFigNames, 1)

%%

% finalize script
clear
close all
