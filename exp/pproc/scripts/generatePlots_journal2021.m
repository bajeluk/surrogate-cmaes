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
supplementaryFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ecml2021journal_supp');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
plotSuppResultsFolder = fullfile(supplementaryFolder, 'images');
tableSuppFolder = fullfile(supplementaryFolder, 'tex');
experimentFolder = fullfile('exp', 'experiments', 'exp_smsp');
resFolder = fullfile(experimentFolder, 'journal');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);
[~, ~] = mkdir(plotSuppResultsFolder);
[~, ~] = mkdir(tableSuppFolder);
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
exp_meta_quantile       = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_quantile.mat'), exp_output, 'Uni', false);
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
mftsNotation = cell(nTSS, 1);
tableMftsNotation = cell(nTSS, 1);
ftsOrdering = cell(nTSS, 1);
ftQuantiles = cell(nTSS, 1);

norm_sigm = @(x, q1, q2, q_val) ...
              1./(1+exp(log((1-q_val)/q_val)* (2*x - q1-q2) / (q1-q2)));

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
  mftsNotation{tss} = ...
    S_sample.mfts_names_red(sort(S_sample.corrMedoidId_nanpair));
  % sort features
  tableMftsNotation{tss} = S_sample.tableMftsNotation_red(sort(S_sample.corrMedoidId_nanpair));
  ftsInside = cellfun(@(x) extractBetween(x, '{', '}'), ...
                      tableMftsNotation{tss}, 'Uni', false);
  [~, ftsOrdering{tss}] = sort(cellfun(@(x) strjoin(x(end:-1:1), ','), ...
                                  ftsInside, 'Uni', false));
  % reorder features
  mftsNotation{tss} = mftsNotation{tss}(ftsOrdering{tss});
  tableMftsNotation{tss} = tableMftsNotation{tss}(ftsOrdering{tss});
  mftsIds = mftsIds(ftsOrdering{tss});

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

  % get quantiles for normalization
  if isfile(exp_meta_quantile{fullId})
    S_q_full = load(exp_meta_quantile{fullId});
    quantileBase = S_q_full.quantile_vals(1);
  else
    error('There is %s missing.', exp_meta_quantile{fullId})
  end
  if (tss == nearId || tss == knnId)
    if isfile(exp_meta_quantile{tss})
      S_q_tss = load(exp_meta_quantile{tss});
    else
      error('There is %s missing.', exp_meta_quantile{tss})
    end
  else
    S_q_tss.mfts_names_red = {};
    S_q_tss.quantiles = [];
  end

  ftIds_full = ismember(S_q_full.mfts_names_red, mftsNotation{tss});
  ftIds_tss  = ismember(S_q_tss.mfts_names_red,  mftsNotation{tss});
  ftQuantiles{tss} = [S_q_full.quantiles(ftIds_full, :); S_q_tss.quantiles(ftIds_tss, :)];
  % fix ordering
  ftQuantiles{tss} = ftQuantiles{tss}(ftsOrdering{tss}, :);
  % normalize metafeatures
  for ft = 1:size(ftQuantiles{tss}, 1)
    % for debugging uncomment the following row
    % fprintf('%42s: %10.2g %10.2g\n', mftsNotation{tss}{ft}, ftQuantiles{tss}(ft, 1), ftQuantiles{tss}(ft, 2))
    mfts{tss}(:, :, :, :, :, ft) = ...
      1./(1+exp(log((1-quantileBase)/quantileBase)* ...
                  (2*mfts{tss}(:, :, :, :, :, ft) - ftQuantiles{tss}(ft, 1)-ftQuantiles{tss}(ft, 2)) / ...
                  (ftQuantiles{tss}(ft, 1)-ftQuantiles{tss}(ft, 2))));
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

% GP and RF settings base original ordering
[gpModelLabelBase, gpModelOrder] = sort({'NN', 'SE', 'LIN', 'Q', ...
                                    'Mat', 'RQ', 'SE+Q', 'Gibbs'});
rfModelLabelBase = {'{CART}_\text{MSE}', '{CART}_\text{RDE}', ...
                 '{SCRT}_\text{MSE}', '{SCRT}_\text{RDE}', ...
                 '{OC1}', ...
                 '{PAIR}_\text{MSE}', '{PAIR}_\text{RDE}', ...
                 '{SUPP}_\text{MSE}', '{SUPP}_\text{RDE}'};
[~, rfModelOrderFull] = sort([1, 2, 6, 4, 3, 8, 7, 5, 9]);
                        % {'Axis_{MSE}', 'Axis_{RDE}', 'Pair_{MSE}', ...
                        % 'Gauss_{RDE}', 'Gauss_{MSE}', 'Res_{MSE}', ...
                        % 'Pair_{RDE}', 'Hill', 'Res_{RDE}'};
[~, rfModelOrderNear] = sort([1, 2, 6, 3, 4, 7, 5, 9, 8]);
                        % {'Axis_{MSE}', 'Axis_{RDE}', 'Pair_{MSE}', ...
                        % 'Gauss_{MSE}', 'Gauss_{RDE}', 'Pair_{RDE}', ...
                        % 'Hill', 'Res_{RDE}', 'Res_{MSE}'};
% add GP and RF model labels
gpModelLabels = cellfun(@(x) ['{{{}}}^\text{', x, '}'], gpModelLabelBase, 'Uni', false);
nGPLabels = numel(gpModelLabels);
rfModelLabels = cellfun(@(x) ['{{}}^\text', x, ''], rfModelLabelBase, 'Uni', false);
nRFLabels = numel(rfModelLabels);
% new labels in sorted order
modelLabels = cell(nTSS, nModels, max(nGPLabels, nRFLabels));
% TSS full GP model
modelLabels(fullId, 1, 1:nGPLabels) = gpModelLabels;
rde(fullId, 1, 1:nGPLabels) = rde(fullId, 1, gpModelOrder);
mse(fullId, 1, 1:nGPLabels) = mse(fullId, 1, gpModelOrder);
% TSS nearest GP model
modelLabels(nearId, 1, 1:nGPLabels) = gpModelLabels;
rde(nearId, 1, 1:nGPLabels) = rde(nearId, 1, gpModelOrder);
mse(nearId, 1, 1:nGPLabels) = mse(nearId, 1, gpModelOrder);
% TSS full RF model
modelLabels(fullId, 2, 1:nRFLabels) = rfModelLabels;
rde(fullId, 1, 1:nRFLabels) = rde(fullId, 1, rfModelOrderFull);
mse(fullId, 1, 1:nRFLabels) = mse(fullId, 1, rfModelOrderFull);
% TSS nearest RF model
modelLabels(nearId, 2, 1:nRFLabels) = rfModelLabels;
rde(nearId, 1, 1:nRFLabels) = rde(nearId, 1, rfModelOrderFull);
mse(nearId, 1, 1:nRFLabels) = mse(nearId, 1, rfModelOrderFull);
% all TSS lmm model
modelLabels(:, 3, 1) = {'{}'};
% all TSS lq model
modelLabels([fullId, nearId], 4, 1) = {'{}'};

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
                  'DataNames', {'GP', 'RF', 'lmm', 'lq', 'GP', 'RF', 'lmm', 'lq', 'lmm'}, ...
                  'DefFile', fullfile(tableFolder, 'defFile.tex'), ...
                  'Evaluations', tableModelLabels, ...
                  'HeadColW', '0.8cm', ...
                  'Mode', 'model', ...
                  'ResultFile', fullfile(tableFolder, ['duelTable_', err{1}, '.tex']), ...
                  'TssList', tssList, ...
                  'TssNums', [19, 19, 1], ...
                  'Vertical', true);
end

%% Table of model failures in prediction
% The number of cases (or percentage) when the model did not provided
% prediction.

percOfNaN = cellfun(@(x) sum(isnan(x(:))), rde)./cellfun(@numel, rde);
modelNums = [8, 9, 1, 1];

% sort model settings
currentLabels = shiftdim(modelLabels(fullId, :, :), 2);
settingsNames = cellfun(@(x) strrep(['$', x, '$'], '\', '\\'), ...
  currentLabels(~cellfun(@isempty, currentLabels(:))), 'Uni', false);

% open file
FID = fopen(fullfile(tableFolder, 'modelNanTable.tex'), 'w');
% print table
fprintf(FID, '\\begin{table}[t]\n');
fprintf(FID, '\n');
% print settings
fprintf(FID, '\\setlength{\\savetabcolsep}{\\tabcolsep}\n');
fprintf(FID, '\\setlength{\\savecmidrulekern}{\\cmidrulekern}\n');
fprintf(FID, '\n');
fprintf(FID, '\\setlength{\\headcolw}{0.85cm}\n');
fprintf(FID, '\\setlength{\\tabcolsep}{0pt}\n');
fprintf(FID, '\\setlength{\\cmidrulekern}{2pt}\n');
fprintf(FID, '\\setlength{\\dueltabcolw}{%0.1f\\textwidth-\\headcolw}\n', ...
             1.0);
fprintf(FID, '\\setlength{\\dueltabcolw}{\\dueltabcolw/%d}\n', sum(modelNums));
fprintf(FID, '\n');
fprintf(FID, '\\centering\n');
% caption
fprintf(FID, '\\caption{\n');
fprintf(FID, ['  Percentages of cases when the model did not provide usable\n', ...
      'prediction (model not trained, its prediction failed, or prediction is constant).']);
fprintf(FID, '}\n\n');
fprintf(FID, '\\label{tab:modelNanTable}\n');
fprintf(FID, '\n');
% table itself
fprintf(FID, '\\begin{tabular}{H%s}\n', repmat('R', 1, sum(modelNums)));
% first top (thick) line
fprintf(FID, '\\toprule\n');
fprintf(FID, ' &');
% header with TSS
modelLine = @(x, y) sprintf('\\\\multicolumn{%d}{l}{\\\\parbox{%d\\\\dueltabcolw}{\\\\centering %s}}', x, x, y{1});
fprintf(FID, strjoin(arrayfun(modelLine, modelNums, {'GP', 'RF', 'lmm', 'lq'}, 'Uni', false), ' & '));
fprintf(FID, '\\\\\n');
% data columns mid-lines
cmidruleCols = [0, cumsum(modelNums), sum(modelNums)];
for i = 2:numel(modelNums)+1
  fprintf(FID, '\\cmidrule(lr){%d-%d}\n', cmidruleCols(i-1)+2, cmidruleCols(i)+1);
end
% header with model settings
fprintf(FID, 'TSS & ');
fprintf(FID, strjoin(settingsNames, ' & ')); % numOfData + 1
fprintf(FID, '\\\\\n');
fprintf(FID, '\\midrule\n');
% percentage values
for tss = 1:nTSS
  fprintf(FID, '%s', tssList{tss});
  for m = 1:numel(modelNums)
    for i = 1:modelNums(m)
      if isnan(percOfNaN(tss, m, i))
        fprintf(FID, ' & ---');
      else
        fprintf(FID, ' & %2.1f', 100*percOfNaN(tss, m, i));
      end
    end
  end
  fprintf(FID, '\\\\\n');
end
fprintf(FID, '\\end{tabular}\n');

fprintf(FID, '\\setlength{\\tabcolsep}{\\savetabcolsep}\n');
fprintf(FID, '\\setlength{\\cmidrulekern}{\\savecmidrulekern}\n');
fprintf(FID, '\\end{table}\n');
% close file
fclose(FID);

clear FID

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

% TSS cycle for TSS full and nearest
% Note: TSS knn included due to reorderedTable variable
for tss = 1:nTSS
  % error cycle
  for err = {'mse', 'rde'}
    ksTabSet.ResultFile = fullfile(tableFolder, ['ksTable_', err{1}, '_', tssList{tss}, '.tex']);
    ksTabSet.ColGroups = {'GP', 'RF', 'lmm', 'lq'};
    ksTabSet.ColGroupNum = [8, 9, 1, 1];
    ksTabSet.HeadColWidth = '2.2cm';
    ksTabSet.GroupHeadColWidth = '0cm';

    currentLabels = modelLabels(tss, :, :);
    
    ksTabSet.Caption = ['The p-values of the Kolmogorov-Smirnov (KS) test comparing the equality of ',...
      'probability distributions of individual TSS ', tssList{tss}, ...
      ' feature representatives on all data ', ...
      'and on those data on which a particular model setting scored best in ', upper(err{1}), ...
      '. ', ...
      'The p-values are after the Holm correction and they are shown only if the KS test ', ...
      'rejects the equality of both distributions ', ...
      'at the family-wise significance level $\alpha=0.05$, ', ...
      'non-rejecting the equality hypothesis is indicated with ---. ', ...
      'Zeros indicate p-values below the smallest double precision number. '];
    ksTabSet.RowValName = ''; % 'TSS';
    
    % get actual table
    if strcmp(err{1}, 'rde')
      actualTable = S_ks.ks_rde_p(:, labelCoor(:, 1) == tss);
    else
      actualTable = S_ks.ks_mse_p(:, labelCoor(:, 1) == tss);
    end
      

    [ksTabSet.ColNames, colNamesId{tss}] = sort(cellfun(@(x) ['$', x, '$'], ...
      currentLabels(~cellfun(@isempty, currentLabels)), 'Uni', false));
    ksTabSet.RowNames = cellfun(@(x) strrep(x, '\', '\\'), ...
                                   tableMftsNotation{tss}, 'Uni', false);
    ksTabSet.RowGroups = {''}; % {tssList{tss}};
    ksTabSet.RowGroupNum = [14];
    ksTabSet.FirstCell = '$\model$';
    ksTabSet.ColValName = 'settings';

    % reorder table
    reorderedTable{tss}.(err{1}) = actualTable(:, colNamesId{tss});
    % print ks table in the right ordering
    prtSignifTable(reorderedTable{tss}.(err{1}), ksTabSet)
    % prepare for image
    mNames{tss} = ksTabSet.ColNames;
  end
end

% table for TSS knn
actualData = mat2cell(...
                [ S_ks.ks_mse_p(:, labelCoor(:, 1) == knnId)';
                  S_ks.ks_rde_p(:, labelCoor(:, 1) == knnId)'], ...
                2, ones(1, 14)...
              );
dataTable = table(actualData{:},...
              'RowNames', {'MSE', 'RDE'}...
            );
lt = LatexTable(dataTable);
lt.setHeaderRow([{''}; cellfun(@(x) ...
        ['\rotatebox{90}{\parbox{2.4cm}{', x, '}}'], ...
        tableMftsNotation{knnId}, 'Uni', false)]);
lt.opts.tableCaption = ...
  ['The p-values of the Kolmogorov-Smirnov (KS) test comparing the equality of ',...
  'probability distributions of individual TSS knn ', ...
  'feature representatives on all data ', ...
  'and on those data on which the lmm model setting scored best in MSE and RDE', ...
  '. ', ...
  'The p-values are after the Holm correction and they are shown only if the KS test ', ...
  'rejects the equality of both distributions ', ...
  'at the family-wise significance level $\alpha=0.05$. ', ...
  'Zeros indicate p-values below the smallest double precision number. '];
lt.opts.booktabs = 1;
lt.opts.tableLabel = 'ksTable_mse_knn';
lt.opts.tableColumnAlignment = repmat({'R'}, 1, 15);
lt.setColumnFormat('\\tiny{%0.2g}')
latexRows = [...
             {'\setlength{\tabcolsep}{0pt}'; ...
              '\newcolumntype{R}{>{\raggedleft\arraybackslash}m{0.75cm}}'}; ...
             lt.toStringRows; ...
             {'\setlength{\tabcolsep}{\savetabcolsep}'}];
% save the result to the file
fid = fopen(fullfile(tableFolder, 'ksTable_mse_knn.tex'), 'w');
for i = 1:length(latexRows)
  fprintf(fid, '%s\n', latexRows{i});
end
fclose(fid);

clear actualTable tss

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
rfModelLabelBase2(strcmp(rfModelLabelBase2, '{OC1}')) = {'{OC1}_\textrm{\hspace{5ex}}'};
gpModelLabels_ks = cellfun(@(x) ['\textrm{', x, '}'], ...
                     sort(gpModelLabelBase), 'Uni', false);
rfModelLabels_ks = cellfun(@(x) ['\textrm', strrep(x, '\text{', '\textrm{'), ''], ...
                     sort(rfModelLabelBase2), 'Uni', false);
simpleModelNames = [gpModelLabels_ks, rfModelLabels_ks, {'\quad', '\quad'}];

% create combinations cov x set
modelSetNames = {};
for m = 1:numel(imgModelLabels)
  for m2 = 1:nModelSettings(m)
    actMName = simpleModelNames{numel(modelSetNames)+1};
    % middle member
    if m2 == round(nModelSettings(m)/2)
      if strcmp(imgModelLabels{m}, '   RF')
        modelSetNames{end+1} = [imgModelLabels{m}, ' \quad$', actMName, '$' ];
      else
        modelSetNames{end+1} = [imgModelLabels{m}, ' \qquad\qquad$', actMName, '$' ];
      end
    % last member gp
    elseif m2 == nModelSettings(m) && m == 1
      modelSetNames{end+1} = ['$\overline{\phantom{AAAA}}\qquad', actMName, '$' ];
    % last member rf
    elseif m2 == nModelSettings(m) && m == 2
      modelSetNames{end+1} = ['$\overline{\phantom{AAAA}}', actMName, '$' ];
    else
      modelSetNames{end+1} = ['$', actMName, '$'];
    end
  end
end

% feature names
imgFeatNames = [];
for tss = 1:nTSS
  imgFeatNames = [imgFeatNames; tableMftsNotation{tss}];
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

%% one plot for TSS full and nearest MSE and RDE

close all

han{1} = figure('Units', 'centimeters', ...
                'Position', [1 1 1.5*sizeX 1.5*sizeY], ...
                'PaperSize', [1.5*sizeX + 2, 1.5*sizeY + 2]);
[rTrows, rTcols] = size(reorderedTable{fullId}.('mse')');
plotData = [reorderedTable{fullId}.('mse')', NaN(rTrows, 1), reorderedTable{nearId}.('mse')'; ...
            NaN(1, 2*rTcols + 1); ...
            reorderedTable{fullId}.('rde')', NaN(rTrows, 1), reorderedTable{nearId}.('rde')'];
          
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
ax.XTickLabel = [imgFeatNames(1:14); {''}; imgFeatNames(15:28)];
ax.YTick = 1:size(plotData, 1);
ax.YTickLabel = [modelSetNames, {''}, modelSetNames];
ax.XTickLabelRotation = labelRot;
ax.TickLabelInterpreter = 'latex';
xlabel('TSS full\hspace{6.5cm}TSS nearest', 'Interpreter', 'latex')
ylabel('RDE\hspace{7.5cm}MSE', 'Interpreter', 'latex')
grid off
% set(ax, 'visible', 'off')
% set(findall(ax, 'type', 'text'), 'visible', 'on')
% line([0,0], [100, 100], 'Color', 'b', 'LineWidth', 10)

% change insignificant data to red
han{1}.Colormap(end, :) = [1 0 0];

hold off

% print figure to pdf
plotNames = {fullfile(plotResultsFolder, 'ks_comb.pdf')};
print2pdf(han, plotNames, 1)

clear ax err han plotNames tss

%% Classification tree analysis
% For TSS full and nearest build a classification tree with models in
% leaves and metafeatures in nodes according to both MSE and RDE. The RDE
% is used as a primary measure. If the value is identical for more than one
% model, MSE is utilized.

fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
         'Performing classification tree analysis\n'], ...
        fix(clock))

% model settings names
settingsNames(1:8) = cellfun(@(x) ['\text{GP}_\text{', x, '}'], gpModelLabelBase, 'Uni', false);
settingsNames(9:17) = cellfun(@(x) ['\text{RF}^\text', x, ''], rfModelLabelBase, 'Uni', false);
settingsNames{18} = '\text{lmm}';
settingsNames{19} = '\text{lq}';

for tss = [fullId, nearId]
  % init
  nMfts = size(mfts{tss}, 6);
  nDataRows = numel(mfts{tss})/nMfts;

  % prepare metafeatures
  mfts_act = reshape(mfts{tss}, [nDataRows, nMfts]);

  % prepare labels
  rde_act = rde(tss, :, :);
  mse_act = mse(tss, :, :);

  existSetId = find(~cellfun(@isempty, rde_act));
  nSettings = numel(existSetId);

  % reshape cell array to double array (last dimension corresponds to
  % individual settings)
  rdeMat = cell2mat(reshape(rde_act(existSetId), [1, 1, 1, 1, 1, nSettings]));
  mseMat = cell2mat(reshape(mse_act(existSetId), [1, 1, 1, 1, 1, nSettings]));
  % reshape again to the number of cases x number of models shape
  rdeMat = reshape(rdeMat, [nDataRows, nSettings]);
  mseMat = reshape(mseMat, [nDataRows, nSettings]);
  % find lowest RDE
  bestSet_rde = min(rdeMat, [], 2);
  % remove rows where lowest RDE is NaN
  mfts_act(isnan(bestSet_rde), :) = [];
  rdeMat(  isnan(bestSet_rde), :) = [];
  mseMat(  isnan(bestSet_rde), :) = [];
  bestSet_rde(isnan(bestSet_rde)) = [];
  % ids of the settings with the lowest RDE
  % get all the mins, where the error is the same
  bestSet_rdeId = rdeMat == repmat(bestSet_rde, [1, nSettings]);
  % replace MSE, where RDE was not lowest with Inf
  mseMat(~bestSet_rdeId) = Inf;
  [~, bestSetId] = min(mseMat, [], 2);

  % low memory limit
  lowMemLimit = min(10^7, numel(bestSetId));
  numOfTreeLevels = 10;
  % tree training
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Tree training for TSS %s start\n'], ...
          fix(clock), tssList{tss})
  CT = fitctree(mfts_act(1:lowMemLimit, :), ...
                settingsNames(bestSetId(1:lowMemLimit)), ...
                'PredictorNames', tableMftsNotation{tss}, ...
                'MinLeafSize', 5000, ...
                'SplitCriterion', 'twoing', ...
                'Categorical', 1);
  % tree pruning
%   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
%            'Pruning trained tree for TSS %s\n'], ...
%           fix(clock), tssList{tss})
%   CT = CT.prune('Level', ...
%             min(max(CT.PruneList), max(CT.PruneList) - numOfTreeLevels));

  % export tree
  settings.LeafProperties = @(a) struct('shape', 'box', 'fillcolor', 'LeafColor', 'style', 'filled');
  settings.LeafFormat = '$%s$';
  settings.InnerNodeProperties = @(a) struct('fillcolor', 'InnerNodeColor', 'style', 'filled');
  settings.InnerNodeFormat = '%s';
  % denormalize thresholds
  settings.CutPoint = CT.CutPoint;
  for cp = 1:numel(CT.CutPredictor)
    ftId = ismember(tableMftsNotation{tss}, CT.CutPredictor(cp));
    if any(ftId)
      settings.CutPoint(cp) = 1/2*(ftQuantiles{tss}(ftId, 1) + ftQuantiles{tss}(ftId, 2) + ...
                                   (ftQuantiles{tss}(ftId, 1) - ftQuantiles{tss}(ftId, 2)) * ...
                                   log(1/CT.CutPoint(cp) - 1) / ...
                                   log((1-quantileBase)/quantileBase));
      settings.CutCategories(cp, :) = cellfun(@(x) ...
                             1/2*(ftQuantiles{tss}(ftId, 1) + ftQuantiles{tss}(ftId, 2) + ...
                                   (ftQuantiles{tss}(ftId, 1) - ftQuantiles{tss}(ftId, 2)) * ...
                                   log(1./x - 1) / ...
                                   log((1-quantileBase)/quantileBase)), CT.CutCategories(cp, :), ...
                                   'Uni', false);
    end
    % observations are natural numbers
    if find(ftId) == 2
      settings.CutPoint(cp) = round(settings.CutPoint(cp));
    else
      % round to 2 numbers after decimal point
      settings.CutPoint(cp) = round(settings.CutPoint(cp)*100)/100;
    end
  end

  % naming edit
  settings.CutPredictor = cellfun(@(x) strrep(x, '\', '\\'), ...
                                      CT.CutPredictor, 'Uni', false);
  treeGVFile = fullfile(plotResultsFolder, ['tree_', tssList{tss}, '_ecml.gv']);
  FID = fopen(treeGVFile, 'w');
  tree2dot(FID, CT, settings);
  fclose(FID);
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Tree for TSS %s saved to %s\n'], ...
          fix(clock), tssList{tss}, treeGVFile)
end

clear cp FID ftId settings tss

%%

% finalize script
clear
close all

fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
         'Script generatePlots_journal2021 finished\n'], ...
         fix(clock))