%% 2021 journal article plots
% Script for analysing surrogate model prediction precision dependence on
% data sampled from DTS-CMA-ES run over the noiseless part of the BBOB
% framework.
% Script also creates graphs and tables.
% 
% Created for 2021 journal article.

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'journal2021_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else

%%

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ecml2021journal');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
experimentFolder = fullfile('exp', 'experiments', 'exp_smsp');
resFolder = fullfile(experimentFolder, 'journal');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);
[~, ~] = mkdir(experimentFolder);
[~, ~] = mkdir(resFolder);

% create tables from calculated features and save them

expfolder = fullfile('exp', 'experiments');
% metafeatures
exp_mfts_id = 'exp_DTSmodels_meta_03_mfts_train';
exp_mfts_knn_id = 'exp_DTSmodels_meta_03_mfts_train_knn';
% models
modelList = {'gp', 'rf', 'lmm', 'lq'};
nModels = numel(modelList);
exp_gp_near_id = 'exp_DTSmodels_meta_03_gp_train';
exp_gp_full_id = 'exp_DTSmodels_meta_03_gp_train_full';
exp_rf_near_id = 'exp_DTSmodels_meta_03_rf_train';
exp_rf_full_id = 'exp_DTSmodels_meta_03_rf_train_full';
exp_lmm_near_id = 'exp_DTSmodels_meta_03_lmm_train';
exp_lmm_full_id = 'exp_DTSmodels_meta_03_lmm_train_full';
exp_lmm_knn_id = 'exp_DTSmodels_meta_03_lmm_train_knn';
exp_lq_near_id = 'exp_DTSmodels_meta_03_lq_train';
exp_lq_full_id = 'exp_DTSmodels_meta_03_lq_train_full';

% TSS settings
tssList = {'full', 'nearest', 'knn'};
nTSS = numel(tssList);
fullId = find(strcmp(tssList, 'full'));
nearId = find(strcmp(tssList, 'nearest'));
knnId  = find(strcmp(tssList, 'knn'));

% results of feature testing
exp_output = cellfun(@(x) fullfile(experimentFolder, x), tssList, 'Uni', false);
exp_smsp_corr_cluster   = cellfun(@(x) fullfile(x, 'exp_smsp_corr_cluster.mat'),   exp_output, 'Uni', false);

% results of statistical testing
exp_journal_friedman = fullfile(resFolder, 'friedman.mat');
exp_journal_kolmogorov = fullfile(resFolder, 'kolmogorov.mat');

exp_mfts_folder = fullfile(expfolder, exp_mfts_id, 'metafeatures');
exp_mfts_knn_folder = fullfile(expfolder, exp_mfts_knn_id, 'metafeatures');

% basic settings
dims = [2, 3, 5, 10, 20];
funs = 1:24;
insts = 11:15;
ids = 1:7;

printScriptMess = false;

clear actualFolder articleFolder

%%

% Load and extract metafeature data

% make list of metafeature files
mftsFileList = searchFile(exp_mfts_folder, '*.mat*');
mftsFileList_knn = searchFile(exp_mfts_knn_folder, '*.mat*');

% init
mfts = cell(nTSS, 1);

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Extracting metafeatures of TSS %s\n'], ...
           fix(clock), tssList{tss})

  % load results of metafeature sample testing (i.e., the ids of
  % metafeatures selected for further research)
  S_sample = load(exp_smsp_corr_cluster{tss});
  % extract ids for metafeatures in sorted order
  mftsIds = sort(S_sample.mftsMedoidIds);
  % extract table names in sorted order
  tableMftsNotation{tss} = ...
    S_sample.mfts_names_red(sort(S_sample.corrMedoidId_nanpair));
  % number of features
  nMfts = numel(mftsIds);
  % extract last archive feature id
  if isfield(S_sample, 'lastArchFtId')
    lastArchFtId = S_sample.lastArchFtId;
  else
    lastArchFtId = 241;
  end
  
  % load metafeatures from experiments on DTS_meta_005_train
  for fl = 1:numel(mftsFileList)
%     fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
%              'Extracting features of TSS %s file %3d/%3d\n'], ...
%              fix(clock), tssList{tss}, fl, numel(mftsFileList))

    % load one file
    S_mf = load(mftsFileList{fl});
    nGen = size(S_mf.values, 2);
    % get fun, dim, inst, and id 
    fId = S_mf.fun == funs;
    dId = S_mf.dim == dims;
    instId = S_mf.inst == insts;
    % id is not saved properly => extract from filename
    idId = str2double(mftsFileList{fl}(strfind(mftsFileList{fl}, '_fts.mat') - 1)) == ids;
    % get values
    vals = S_mf.values;
    if tss == knnId
      S_mf_knn = load(mftsFileList_knn{fl});
      vals(lastArchFtId+1:end, :) = S_mf_knn.values;
    end
    % cat file results to the rest
    mfts{tss}(dId, fId, instId, idId, 1:nGen, 1:nMfts) = ...
      reshape(vals(mftsIds, :)', [1, 1, 1, 1, nGen, nMfts]);
  end
end

clear('mftsFileList', 'mftsFileList_knn', 'tss S_sample', 'mftsIds', ...
      'nMfts', 'lastArchFtId', 'fl', 'S_mf', 'fId', 'dId', 'instId', ...
      'idId', 'vals')

%%

% Load calculated data

if printScriptMess
  fprintf('Loading results\n')
end

% init
rde = cell(nTSS, nModels);
mse = cell(nTSS, nModels);

% TSS method loop
for tss = 1:nTSS
  % TSS method addition
  if tss == nearId
    addTSS = [];
  else
    addTSS = ['_', tssList{tss}];
  end
  % model type loop
  for mt = 1:nModels
    modelFolder = fullfile(expfolder, ['exp_DTSmodels_meta_03_', ...
      modelList{mt}, '_train', addTSS]);
    if isdir(modelFolder)
      % find different model settings in the model folder
      settingsDir = dir(fullfile(modelFolder, '*model_*'));
      settingsList = {settingsDir.name};
      nModelSettings = numel(settingsList);
      % model settings loop
      for ms = 1:nModelSettings
        fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
          'Loading results of TSS %s model %s settings %d\n'], ...
          fix(clock), tssList{tss}, modelList{mt}, ms)

        fileList = searchFile(fullfile(modelFolder, settingsList{ms}), '*.mat');
        nFiles = numel(fileList);
        % init error values
        rde{tss, mt, ms} = NaN(numel(dims), numel(funs), numel(insts), numel(ids), nGen);
        mse{tss, mt, ms} = NaN(numel(dims), numel(funs), numel(insts), numel(ids), nGen);
        % file loop
        for fl = 1:nFiles
          % extract results
          S = load(fileList{fl});
          if fl == 1
            switch modelList{mt}
              case 'gp'
                gpSetOrder{tss, ms} = S.modelOptions1.hypOptions.covFcn;
              case 'rf'
                rfSetOrder{tss, ms} = S.modelOptions1;
              case 'lmm'
                lmmSetOrder{tss, ms} = S.modelOptions1;
              case 'lq'
                lqSetOrder{tss, ms} = S.modelOptions1;
            end
          end
          % get fun, dim, inst, and id
          fId = S.fun == funs;
          dId = S.dim == dims;
          instId = S.instances == insts;
          idId = S.ids == ids;
          % get number of generations
          nGen = size(S.stats.rde, 3);
          % get RDE values
          rde{tss, mt, ms}(dId, fId, instId, idId, 1:nGen) = ...
            reshape(S.stats.rde, [1, 1, sum(instId), sum(idId), nGen]);
          % get MSE values
          mse{tss, mt, ms}(dId, fId, instId, idId, 1:nGen) = ...
            reshape(S.stats.mse, [1, 1, sum(instId), sum(idId), nGen]);
        end
      end
    end
  end
end

clear('addTSS', 'modelFolder', 'mt', 'settingsDir', 'settingsList', 'tss', ...
      'nModelSettings', 'ms', 'fileList', 'nFiles', 'fl', 'S', 'fId', ...
      'dId', 'instId', 'idId')

%%

% create model labels

modelLabels = cell(size(rde));
% base GP and RF settings labels
gpModelLabelBase = {'NN', 'SE', 'LIN', 'Q', ...
                 'Mat', 'RQ', 'SE+Q', 'Gibbs'};
rfModelLabelBase = {'{Axis}_\text{MSE}', '{Axis}_\text{RDE}', ...
                 '{Gauss}_\text{MSE}', '{Gauss}_\text{RDE}', ...
                 '{Hill}', ...
                 '{Pair}_\text{MSE}', '{Pair}_\text{RDE}', ...
                 '{Res}_\text{MSE}', '{Res}_\text{RDE}'};
% add GP and RF model labels
gpModelLabels = cellfun(@(x) ['{{{}}}^\text{', x, '}'], gpModelLabelBase, 'Uni', false);
rfModelLabels = cellfun(@(x) ['{{}}^\text', x, ''], rfModelLabelBase, 'Uni', false);
% TSS full GP model
modelLabels(1, 1, 1:8) = gpModelLabels;
% TSS nearest GP model
modelLabels(2, 1, 1:8) = gpModelLabels;
% TSS full RF model
modelLabels(1, 2, 1:9) = rfModelLabels([1, 2, 6, 4, 3, 8, 7, 5, 9]);
                        % {'Axis_{MSE}', 'Axis_{RDE}', 'Pair_{MSE}', ...
                        % 'Gauss_{RDE}', 'Gauss_{MSE}', 'Res_{MSE}', ...
                        % 'Pair_{RDE}', 'Hill', 'Res_{RDE}'};
% TSS nearest RF model
modelLabels(2, 2, 1:9) = rfModelLabels([1, 2, 6, 3, 4, 7, 5, 9, 8]);
                        % {'Axis_{MSE}', 'Axis_{RDE}', 'Pair_{MSE}', ...
                        % 'Gauss_{MSE}', 'Gauss_{RDE}', 'Pair_{RDE}', ...
                        % 'Hill', 'Res_{RDE}', 'Res_{MSE}'};
% all TSS lmm model
modelLabels(1:3, 3, 1) = {'{}'};
% all TSS lq model
modelLabels(1:2, 4, 1) = {'{}'};

% add TSS labels
for tss = 1:nTSS
  existLabel = ~cellfun(@isempty, modelLabels(tss, :, :));
%   modelLabels(tss, existLabel) = cellfun(@(x) [x, '^{', tssList{tss}, '}'], ...
%                                      modelLabels(tss, existLabel), 'Uni', false);
end

clear('tss')

%%
if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Test pairwise differences in models’ errors
% Rejecting independence of models' errors using Friedman test with
% Bonferroni-Holm correction on the alpha level.

alpha = 0.05;

% init
existSetId = find(~cellfun(@isempty, modelLabels));
settingsLabels = modelLabels(existSetId);
nSettings = numel(existSetId);
settingCombs = cell(nchoosek(nSettings, 2), 2);
wilcoxon_rde_p = NaN(nSettings);
wilcoxon_mse_p = NaN(nSettings);
wilcoxon_rde_p_bh = NaN(nSettings);
wilcoxon_mse_p_bh = NaN(nSettings);
duel_rde = NaN(nSettings);
duel_mse = NaN(nSettings);

fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
         'Performing Friedman test on errors\n'], ...
        fix(clock))

% prepare data
rde_data = [];
mse_data = [];
% cat errors for individual model settings
for ms = 1:nSettings
  rde_data = [rde_data, rde{existSetId(ms)}(:)];
  mse_data = [mse_data, mse{existSetId(ms)}(:)];
end

% remove NaN's
rde_data(any(isnan(rde_data), 2), :) = [];
mse_data(any(isnan(mse_data), 2), :) = [];
if size(rde_data, 1) > 1 && size(rde_data, 2) > 1
  % Friedman's test on statistical significance of pairwise differences
  % of RDE
  [friedman_rde_p, ~, friedman_stats_rde] = friedman(rde_data, 1, 'off');
  % Friedman's test on statistical significance of pairwise differences
  % of MSE
  [friedman_mse_p, ~, friedman_stats_mse] = friedman(mse_data, 1, 'off');
  % calculate Wilcoxon on pairs when rejected
  if friedman_rde_p < alpha || friedman_mse_p < alpha
    % Wilcoxon's cycle
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
             'Performing two-sided Wilcoxon signed rank test on pairs of model settings errors\n'], ...
            fix(clock))
    j = 0;
    for ms = 1:nSettings
      for ms2 = ms+1:nSettings
        % increase counter
        j = j+1;
        % prepare pair data
        rde_pair = [rde{existSetId(ms)}(:), rde{existSetId(ms2)}(:)];
        mse_pair = [mse{existSetId(ms)}(:), mse{existSetId(ms2)}(:)];
        pairRowIds_rde = all(~isnan(rde_pair), 2);
        pairRowIds_mse = all(~isnan(mse_pair), 2);
        rde_1 = rde_pair(pairRowIds_rde, 1);
        rde_2 = rde_pair(pairRowIds_rde, 2);
        mse_1 = mse_pair(pairRowIds_mse, 1);
        mse_2 = mse_pair(pairRowIds_mse, 2);
        % Wilcoxon test on all non-nan error values of a pair of settings
        wilcoxon_rde_p(ms, ms2) = signrank(rde_1, rde_2);
        wilcoxon_mse_p(ms, ms2) = signrank(mse_1, mse_2);
        % number of wins in lower errors
        duel_rde(ms, ms2) = sum(rde_1 < rde_2);
        duel_rde(ms2, ms) = sum(rde_1 > rde_2);
        duel_mse(ms, ms2) = sum(mse_1 < mse_2);
        duel_mse(ms2, ms) = sum(mse_1 > mse_2);
        % save settings combinations
        settingCombs(j, :) = {modelLabels{existSetId(ms)}, modelLabels{existSetId(ms2)}};
      end
    end
    % p-values using Bonferroni-Holm correction on the alpha level
    [wilcoxon_rde_p_bh(~isnan(wilcoxon_rde_p)), wilcoxon_rde_bh_h] = ...
      bonfHolm(wilcoxon_rde_p(~isnan(wilcoxon_rde_p)), alpha);
    [wilcoxon_mse_p_bh(~isnan(wilcoxon_mse_p)), wilcoxon_mse_bh_h] = ...
      bonfHolm(wilcoxon_mse_p(~isnan(wilcoxon_mse_p)), alpha);
  end
end

% not significant pair results
acceptedPairs_rde = settingCombs(~wilcoxon_rde_bh_h, :);
acceptedPairs_mse = settingCombs(~wilcoxon_mse_bh_h, :);
% make matrices with p-values symmetrical
wilcoxon_rde_p = triu(wilcoxon_rde_p) + triu(wilcoxon_rde_p)';
wilcoxon_mse_p = triu(wilcoxon_mse_p) + triu(wilcoxon_mse_p)';
wilcoxon_rde_p_bh = triu(wilcoxon_rde_p_bh) + triu(wilcoxon_rde_p_bh)';
wilcoxon_mse_p_bh = triu(wilcoxon_mse_p_bh) + triu(wilcoxon_mse_p_bh)';
% calculate percentages of duels
duel_rde_perc = duel_rde./(duel_rde + duel_rde');
duel_mse_perc = duel_mse./(duel_mse + duel_mse');
  
% save testing results
save(exp_journal_friedman, ...
    'alpha', 'settingsLabels', 'settingCombs', ...
    'friedman_rde_p', 'friedman_stats_rde', ...
    'friedman_mse_p', 'friedman_stats_mse', ...
    'wilcoxon_rde_p', 'wilcoxon_rde_p_bh', 'wilcoxon_rde_bh_h', ...
    'wilcoxon_mse_p', 'wilcoxon_mse_p_bh', 'wilcoxon_mse_bh_h', ...
    'acceptedPairs_rde', 'acceptedPairs_mse', ...
    'duel_rde', 'duel_mse', 'duel_rde_perc', 'duel_mse_perc')

% clear variables saved in exp_journal_friedman
clear('alpha', 'j', 'ms', 'mse_data', 'rde_data', ...
    'friedman_rde_p', 'friedman_stats_rde', ...
    'friedman_mse_p', 'friedman_stats_mse', ...
    'wilcoxon_rde_p', 'wilcoxon_rde_p_bh', 'wilcoxon_rde_bh_h', ...
    'wilcoxon_mse_p', 'wilcoxon_mse_p_bh', 'wilcoxon_mse_bh_h', ...
    'rde_pair', 'mse_pair', 'pairRowIds_rde', 'pairRowIds_mse', ...
    'rde_1', 'rde_2', 'mse_1', 'mse_2', ...
    'acceptedPairs_rde', 'acceptedPairs_mse', ...
    'duel_rde', 'duel_mse', 'duel_rde_perc', 'duel_mse_perc')

%% Table of pairwise differences in model error's

% load WCX test results
if isfile(exp_journal_friedman)
  S_wcx = load(exp_journal_friedman);
else
  error('File %s is missing!', exp_journal_friedman)
end

% ordering according to TSS
labelTSSOrder = [1:2:6, 8:2:39, ...
                 2:2:6, 9:2:39, ...
                 7];
tableModelLabelsBase = S_wcx.settingsLabels(labelTSSOrder);

tableModelLabels = {};
tableModelLabelsId = [];
for tss = 1:nTSS
  labCoor = ((tss-1)*19+1) : min(tss*19, numel(labelTSSOrder));
  currentLabels = tableModelLabelsBase(labCoor);
  [  tableModelLabels(labCoor), tableModelLabelsId(labCoor)] = ...
     sort(cellfun(@(x) ['$', x, '$'], currentLabels, 'Uni', false));
   tableModelLabelsId(labCoor) = tableModelLabelsId(labCoor) + 19*(tss-1);
end
% prepare for sprintf
tableModelLabels = cellfun(@(x) strrep(x, '\', '\\'), ...
                     tableModelLabels, 'Uni', false);

% error cycle
for err = {'mse', 'rde'}
  tableData = round(100*S_wcx.(['duel_', err{1}, '_perc']));
  % add decimals to 50% values
  tableData(tableData == 50) = 100*S_wcx.(['duel_', err{1}, '_perc'])(tableData == 50);
  % order table data
  tableData = tableData(labelTSSOrder(tableModelLabelsId), labelTSSOrder(tableModelLabelsId));
  pvals = S_wcx.(['wilcoxon_', err{1}, '_p_bh'])(labelTSSOrder(tableModelLabelsId), labelTSSOrder(tableModelLabelsId));

  % print duel table
  wcxDuelTable(tableData, pvals, ...
                  'DataCaption', upper(err{1}), ...
                  'DataDims', 2, ...
                  'DataNames', {'GP', 'RF', 'lq', 'lmm', 'GP', 'RF', 'lq', 'lmm', 'lmm'}, ...
                  'DefFile', fullfile(tableFolder, 'defFile.tex'), ...
                  'Evaluations', tableModelLabels, ...
                  'HeadColW', '0.8cm', ...
                  'Mode', 'model', ...
                  'ResultFile', fullfile(tableFolder, ['duelTable_', err{1}, '.tex']), ...
                  'TssList', tssList, ...
                  'TssNums', [19, 19, 1], ...
                  'Vertical', true);
end

%% Kolmogorov-Smirnov test
% Kolmogorov-Smirnov test testing equality of distribution of an individual
% feature values calculated on the datasets where the model error is the
% lowest for setting among all tested settings and their distribution on
% all datasets.

% init
existSetId = find(~cellfun(@isempty, modelLabels));
nSettings = numel(existSetId);

% reshape cell array to double array (last dimension corresponds to
% individual settings)
rdeMat = cell2mat(reshape(rde(existSetId), [1, 1, 1, 1, 1, nSettings]));
mseMat = cell2mat(reshape(mse(existSetId), [1, 1, 1, 1, 1, nSettings]));
% ids of the settings with the lowest RDE
% get all the mins, where the error is the same
bestSet_rde = min(rdeMat, [], 6);
bestSet_mse = min(mseMat, [], 6);

existSetCoor = ind2coor(existSetId, size(modelLabels));

bestSetNum_rde = NaN(1, nSettings);
bestSetNum_mse = NaN(1, nSettings);
ks_rde_p = NaN(size(mfts{1}, 6), nSettings);
ks_mse_p = NaN(size(mfts{1}, 6), nSettings);
% KS testing RDE
for s = 1:nSettings
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing Kolmogorov-Smirnov test on metafeatures for model %2d/%d\n'], ...
          fix(clock), s, nSettings)
  tssId = existSetCoor(s, 1, 1);
  bestSetId_rde = rdeMat(:, :, :, :, :, s) == bestSet_rde;
  bestSetId_mse = mseMat(:, :, :, :, :, s) == bestSet_mse;
  bestSetNum_rde(s) = numel(find(bestSetId_rde));
  bestSetNum_mse(s) = numel(find(bestSetId_mse));
  % feature loop
  for m = 1:size(mfts{1}, 6)
    mftsData = mfts{tssId}(:, :, :, :, :, m);
    [~, ks_rde_p(m, s)] = kstest2(mftsData(bestSetId_rde), mftsData(:));
    [~, ks_mse_p(m, s)] = kstest2(mftsData(bestSetId_mse), mftsData(:));
  end
end

save(exp_journal_kolmogorov, 'ks_rde_p', 'ks_mse_p', ...
                             'bestSetNum_rde', 'bestSetNum_mse')

clear('rdeMat', 'mseMat', ...
      'bestSet_rde', 'bestSet_mse', ...
      'bestSetNum_rde', 'bestSetNum_mse', ...
      'bestSetId_rde', 'bestSetId_mse', ...
      'm', 'mftsData''s', 'tssId')

%%
% KS test table

% load KS test results
if isfile(exp_journal_kolmogorov)
  S_ks = load(exp_journal_kolmogorov, 'ks_rde_p', 'ks_mse_p');
else
  error('File %s is missing!', exp_journal_kolmogorov)
end

% prepare label coordinates
labelCoor = ind2coor(find(~cellfun(@isempty, modelLabels)), size(modelLabels));

  % TSS cycle
for tss = 1:nTSS
  if isfile(exp_smsp_corr_cluster{tss})
    S_cl = load(exp_smsp_corr_cluster{tss}, ...
                'tableMftsNotation', 'corrMedoidId_nanpair');
  else
    error('File %s is missing!', exp_smsp_corr_cluster{tss})
  end

  % error cycle
  for err = {'mse', 'rde'}
    ksTabSet.ResultFile = fullfile(tableFolder, ['ksTable_', err{1}, '_', tssList{tss}, '.tex']);
    if tss == knnId
      ksTabSet.ColGroups = {'lmm'};
      ksTabSet.ColGroupNum = 1;
    else
      ksTabSet.ColGroups = {'GP', 'RF', 'lmm', 'lq'};
      ksTabSet.ColGroupNum = [8, 9, 1, 1];
    end
    ksTabSet.ColValName = 'settings';
    ksTabSet.HeadColWidth = '2.4cm';

    currentLabels = modelLabels(tss, :, :);
    [ksTabSet.ColNames, colNamesId{tss}] = sort(cellfun(@(x) ['$', x, '$'], ...
      currentLabels(~cellfun(@isempty, currentLabels)), 'Uni', false));
    % sort features
    tabMftsNot{tss} = S_cl.tableMftsNotation(S_cl.corrMedoidId_nanpair);
    ftsInside = cellfun(@(x) extractBetween(x, '{', '}'), ...
                        tabMftsNot{tss}, 'Uni', false);
    [~, ftsOrdering{tss}] = sort(cellfun(@(x) strjoin(x(end:-1:1), ','), ...
                                    ftsInside, 'Uni', false));
    ksTabSet.RowNames = cellfun(@(x) strrep(x, '\', '\\'), ...
                                   tabMftsNot{tss}(ftsOrdering{tss}), 'Uni', false);
    % TODO: group features according to classes
    ksTabSet.RowGroups = {tssList{tss}};
    ksTabSet.RowGroupNum = [14];
    ksTabSet.RowValName = 'TSS';
    ksTabSet.FirstCell = '$\model$';

    ksTabSet.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality of ',...
          'probability distributions of individual TSS ', tssList{tss}, ...
          ' feature representatives on all data ', ...
          'and on those data on which a particular model setting scored best in ', upper(err{1}), ...
          '. ', ...
          'P-values denote KS test results ', ...
          'rejecting the equality of both distributions with the Holm ', ...
          'correction at the family-wise significance level $\alpha=0.05$, ', ...
          'otherwise, p-values are not shown. ', ...
          'Zeros indicate p-values below the smallest double precision number. '];%, ...
          % '$\featMM$ notation: \texttt{l} -- \texttt{lin}, \texttt{q} -- \texttt{quad}, \texttt{s} -- \texttt{simple}, ', ...
          % '\texttt{i} -- \texttt{interact}, \texttt{c} -- \texttt{coef}.'];

    % get actual table
    if strcmp(err{1}, 'rde')
      actualTable = S_ks.ks_rde_p(:, labelCoor(:, 1) == tss);
    else
      actualTable = S_ks.ks_mse_p(:, labelCoor(:, 1) == tss);
    end
    % reorder table
    reorderedTable{tss}.(err{1}) = actualTable(ftsOrdering{tss}, colNamesId{tss});
    % print ks table in the right ordering
    prtSignifTable(reorderedTable{tss}.(err{1}), ksTabSet)
    % prepare for image
    mNames{tss} = ksTabSet.ColNames;
  end
end

%% 
% KS test image
close all

% image settings
sizeY = 16;
sizeX = 16;
sizeYknn = 8;
sizeXknn = 6.7;
labelRot = 60;

% image model labels
imgModelLabels = {'   GP', '   RF', '  lmm', '    lq'};
nModelSettings = [8, 9, 1, 1];
rfModelLabelBase2 = rfModelLabelBase;
rfModelLabelBase2(strcmp(rfModelLabelBase2, 'Hill')) = {'{Hill}_\textrm{\ \ \ }'};
gpModelLabels_ks = cellfun(@(x) ['\textrm{', x, '}'], ...
                     sort(gpModelLabelBase), 'Uni', false);
rfModelLabels_ks = cellfun(@(x) ['\textrm', strrep(x, '\text{', '\textrm{'), ''], ...
                     sort(rfModelLabelBase), 'Uni', false);
simpleModelNames = [gpModelLabels_ks, rfModelLabels_ks, {'\quad', '\quad'}];

% create combinations cov x set
modelSetNames = {};
for m = 1:numel(imgModelLabels)
  for m2 = 1:nModelSettings(m)
    actMName = simpleModelNames{numel(modelSetNames)+1};
    % middle member
    if m2 == round(nModelSettings(m)/2)
      modelSetNames{end+1} = [imgModelLabels{m}, ' \qquad\qquad$', actMName, '$' ];
    % last member gp
    elseif m2 == nModelSettings(m) && m == 1
      modelSetNames{end+1} = ['$\overline{\phantom{AAAA}}\qquad', actMName, '$' ];
    % last member rf
    elseif m2 == nModelSettings(m) && m == 2
      modelSetNames{end+1} = ['$\overline{\phantom{AAAA}}\quad', actMName, '$' ];
    else
      modelSetNames{end+1} = ['$', actMName, '$'];
    end
  end
end

% feature names
imgFeatNames = [];
for tss = 1:nTSS
  imgFeatNames = [imgFeatNames; tabMftsNot{tss}(ftsOrdering{tss})];
end
% create labels identical to the rest of the paper
imgFeatNamesParts = extractBetween(imgFeatNames, '{', '}');
for i = 1:numel(imgFeatNames)
  imgFeatNames{i} = ['$\varphi_\texttt{', imgFeatNamesParts{i, 1}, '}'];
  switch imgFeatNamesParts{i, 2}
    case '\archive'
      imgFeatNames{i} = [imgFeatNames{i}, '(\mathcal{A}'];
    case '\archivepred'
      imgFeatNames{i} = [imgFeatNames{i}, '(\mathcal{A}_\mathcal{P}'];
    case '\trainset'
      imgFeatNames{i} = [imgFeatNames{i}, '(\mathcal{T}'];
    case '\trainpredset'
      imgFeatNames{i} = [imgFeatNames{i}, '(\mathcal{T}_\mathcal{P}'];
  end
  if strcmp(imgFeatNamesParts{i, 3}, '\transcma')
    imgFeatNames{i} = [imgFeatNames{i}, '^{{}_\top}'];
  end
  if ~isempty(imgFeatNamesParts{i, 2})
    imgFeatNames{i} = [imgFeatNames{i}, ')'];
  end
  imgFeatNames{i} = [imgFeatNames{i}, '$'];
end
% create same lengths
% imgFeatNames = cellfun(@(x) sprintf('%30s', x), imgFeatNames, 'UniformOutput', false);

for err = {'mse', 'rde'}
  for tss = 1:nTSS
    % draw coeffs without colorbar
    if tss == knnId
      han{1} = figure('Units', 'centimeters', ...
                      'Position', [1 1 sizeXknn sizeYknn], ...
                      'PaperSize', [sizeXknn + 2, sizeYknn + 2]);
    else
      han{1} = figure('Units', 'centimeters', ...
                      'Position', [1 1 sizeX sizeY], ...
                      'PaperSize', [sizeX + 2, sizeY + 2]);
    end

    % adjust data
    if tss == knnId
      plotData = reorderedTable{tss}.(err{1});
    else
      plotData = reorderedTable{tss}.(err{1})';
    end

    % non-significant data
    inSignifId = plotData > 0.05;
    % data closest to significance level
    secondSignifVal = 1e-6;
    closeData = plotData > secondSignifVal & plotData <= 0.05;
    % change values to preserve one color for insignificant data
    plotData(closeData) = secondSignifVal;

    % draw image
    hold on
    imagesc(log10(plotData), 'AlphaData',~isnan(plotData))
    colorbar('Ticks', [-300, -250, -200, -150, -100, -50, -5],...
             'TickLabels', {'10^{-300}', '10^{-250}', '10^{-200}', ...
                            '10^{-150}', '10^{-100}', '10^{-50}', '0.05'})

    % axis square
    ax = gca;
    ax.XTick = 1:size(plotData, 2);
    if tss == knnId
      ax.XTickLabel = '';
    else
      ax.XTickLabel = imgFeatNames((tss-1)*14 + (1:14));
    end
    ax.YTick = 1:size(plotData, 1);
    if tss == knnId
      ax.YTickLabel = imgFeatNames((tss-1)*14 + (1:14));
    else
      ax.YTickLabel = modelSetNames;
    end
    ax.XTickLabelRotation = labelRot;
    ax.TickLabelInterpreter = 'latex';
    if tss == knnId
      xlabel(upper(err{1}), 'Interpreter', 'latex')
      ylabel('TSS knn', 'Interpreter', 'latex')
    else
      ylabel(upper(err{1}), 'Interpreter', 'latex')
    end
    grid off
    % set(ax, 'visible', 'off')
    % set(findall(ax, 'type', 'text'), 'visible', 'on')

    % change insignificant data to red
    han{1}.Colormap(end, :) = [1 0 0];

    hold off

    % print figure to pdf
    plotNames = {fullfile(plotResultsFolder, ...
                   sprintf('ks_%s_%s.pdf', tssList{tss}, err{1}))};
    print2pdf(han, plotNames, 1)
  end
end

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
        ylabel('PDF', 'Interpreter', 'latex')
      end
%       title([ks_res.MetafeatureNames{r, c}, ' for ', modelLabels{m}])
      title(modelLabels{m}, 'Interpreter', 'latex')
      legend({'all', modelLabels{m}}, 'Interpreter', 'latex')
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
