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
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'smsppaper');
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
gpModelLabels = {'NN', 'SE', 'LIN', 'Q', ...
                 'Mat', 'RQ', 'SE+Q', 'Gibbs'};
rfModelLabels = {'Axis_{MSE}', 'Axis_{RDE}', ...
                 'Gauss_{MSE}', 'Gauss_{RDE}', ...
                 'Hill', ...
                 'Pair_{MSE}', 'Pair_{RDE}', ...
                 'Res_{MSE}', 'Res_{RDE}'};
% add GP and RF model labels
gpModelLabels = cellfun(@(x) ['gp_', x], gpModelLabels, 'Uni', false);
rfModelLabels = cellfun(@(x) ['rf_', x], rfModelLabels, 'Uni', false);
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
modelLabels(1:3, 3, 1) = {'lmm'};
% all TSS lq model
modelLabels(1:2, 4, 1) = {'lq'};

% add TSS labels
for tss = 1:nTSS
  existLabel = ~cellfun(@isempty, modelLabels(tss, :, :));
  modelLabels(tss, existLabel) = cellfun(@(x) [x, '^{', tssList{tss}, '}'], ...
                                     modelLabels(tss, existLabel), 'Uni', false);
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
wilcoxon_rde_p = NaN(1, nchoosek(nSettings, 2));
wilcoxon_mse_p = NaN(1, nchoosek(nSettings, 2));

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
        rde_pair = [rde{existSetId(ms)}(:), rde{existSetId(ms2)}(:)];
        mse_pair = [mse{existSetId(ms)}(:), mse{existSetId(ms2)}(:)];
        % Wilcoxon test on all non-nan error values of a pair of settings
        wilcoxon_rde_p(1, j) = signrank(rde_pair(all(~isnan(rde_pair), 2), 1), ...
                                        rde_pair(all(~isnan(rde_pair), 2), 2));
        wilcoxon_mse_p(1, j) = signrank(mse_pair(all(~isnan(mse_pair), 2), 1), ...
                                        mse_pair(all(~isnan(mse_pair), 2), 2));
        % save settings combinations
        settingCombs(j, :) = {modelLabels{existSetId(ms)}, modelLabels{existSetId(ms2)}};
      end
    end
    % p-values using Bonferroni-Holm correction on the alpha level
    wilcoxon_rde_p_bh = bonfHolm(wilcoxon_rde_p, alpha);
    wilcoxon_mse_p_bh = bonfHolm(wilcoxon_rde_p, alpha);
  end
end

acceptedPairs_rde = settingCombs(S.wilcoxon_rde_p_bh > alpha, :);
acceptedPairs_mse = settingCombs(S.wilcoxon_mse_p_bh > alpha, :);
  
% save testing results
save(exp_journal_friedman, ...
    'alpha', 'settingsLabels', 'settingCombs', ...
    'friedman_rde_p', 'friedman_stats_rde', 'wilcoxon_rde_p', 'wilcoxon_rde_p_bh', ...
    'friedman_mse_p', 'friedman_stats_mse', 'wilcoxon_mse_p', 'wilcoxon_mse_p_bh', ...
    'acceptedPairs_rde', 'acceptedPairs_mse')

% clear variables saved in exp_journal_friedman
clear('alpha', 'j', 'ms', 'mse_data', 'rde_data', ...
    'friedman_rde_p', 'friedman_stats_rde', 'wilcoxon_rde_p', 'wilcoxon_rde_p_bh', ...
    'friedman_mse_p', 'friedman_stats_mse', 'wilcoxon_mse_p', 'wilcoxon_mse_p_bh', ...
    'rde_pair', 'mse_pair', ...
    'acceptedPairs_rde', 'acceptedPairs_mse')
  
%% Kolmogorov-Smirnov test
% Kolmogorov-Smirnov test testing equality of distribution of an individual
% feature values calculated on the datasets where the model error is the
% lowest for setting among all tested settings and their distribution on
% all datasets.

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

  ksRdeTabSet.ResultFile = fullfile(tableFolder, ['ksTable_rde_', tssList{tss}, '.tex']);
  ksMseTabSet.ResultFile = fullfile(tableFolder, ['ksTable_mse_', tssList{tss}, '.tex']);
  ksRdeTabSet.ColGroups = {'GP', 'RF', 'lmm', 'lq'};
  ksRdeTabSet.ColGroupNum = [8, 9, 1, 1];
  ksRdeTabSet.ColValName = 'settings';
  % tssLabelIds = find 
  currentLabels = modelLabels(tss, :, :);
  ksRdeTabSet.ColNames = currentLabels(~cellfun(@isempty, currentLabels));
  ksRdeTabSet.RowNames = S_cl.tableMftsNotation(S_cl.corrMedoidId_nanpair);
  % TODO: group features according to classes
  ksRdeTabSet.RowGroups = {'fts'};
  ksRdeTabSet.RowGroupNum = [42];
  ksRdeTabSet.RowValName = 'fts group';
 
  ksRdeTabSet.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality of ',...
        'probability distributions of individual features on all data ', ...
        'and on those data on which a particular covariance function scored best. ', ...
        'P-values denote KS test results ', ...
        'rejecting the equality of both distributions with the Holm ', ...
        'correction at the family-wise significance level $\alpha=0.05$, ', ...
        'otherwise, p-values are not shown. ', ...
        '''{\normalfont X}'' values in $\trainpredset$ column denote features impossible to calculate with points without fitness values. ', ...
        'Zeros indicate p-values below the smallest double precision number. ', ...
        '$\featMM$ notation: \texttt{l} -- \texttt{lin}, \texttt{q} -- \texttt{quad}, \texttt{s} -- \texttt{simple}, ', ...
        '\texttt{i} -- \texttt{interact}, \texttt{c} -- \texttt{coef}.'];

  % print ks table
  prtSignifTable(S_ks.ks_rde_p(:, labelCoor(:, 1) == tss), ksRdeTabSet)
end

%%
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
ksTabSet.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality of ',...
        'probability distributions of individual features on all data ', ...
        'and on those data on which a particular covariance function scored best. ', ...
        'P-values denote KS test results ', ...
        'rejecting the equality of both distributions with the Holm ', ...
        'correction at the family-wise significance level $\alpha=0.05$, ', ...
        'otherwise, p-values are not shown. ', ...
        '''{\normalfont X}'' values in $\trainpredset$ column denote features impossible to calculate with points without fitness values. ', ...
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
ksTabSetInd.Caption = ['The p-values of the Kolmogorov-Smirnov test comparing the equality of ',...
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

% image model labels
imgModelLabels = {'   NN', '   SE', '  LIN', '    Q', '  Mat', '   RQ', ' SE+Q', 'Gibbs'};
% set names
setNames = {'$\mathcal{A}\ $', '$\mathcal{T}\ $', '\mathcal{T}_\mathcal{P}'};
% create combinations cov x set
modelSetNames = {};
for m = 1:nModel
  modelSetNames{end+1} = setNames{1};
  modelSetNames{end+1} = [imgModelLabels{m}, ' \ \ ', setNames{2} ];
  if m < nModel
    modelSetNames{end+1} = ['$\overline{\phantom{AAA}}\quad', setNames{3}, '$' ];
  else
    modelSetNames{end+1} = ['$' setNames{3} '$'];
  end
%   end
end

% feature names
imgFeatNames = { ...
    '$N$ \quad\texttt{observations}', ...
    ... cmaes
    '$\overline{\phantom{,}\Phi_\textrm{CMA}}$  \texttt{cma\_mean\_dist}', ...
    '\texttt{cma\_lik}', ...
    ... dispersion
    '$\overline{\phantom{AAAA}}$\quad\texttt{ratio\_mean\_02}', ...
    '\texttt{ratio\_median\_02}', ...
    '\texttt{diff\_mean\_02}', ...
    '\texttt{diff\_median\_02}', ...
    '\texttt{ratio\_mean\_05}', ...
    '\texttt{ratio\_median\_05}', ...
    '\texttt{diff\_mean\_05}', ...
    '$\Phi_\textrm{Dis}$\ \ \texttt{diff\_median\_05}', ...
    '\texttt{ratio\_mean\_10}', ...
    '\texttt{ratio\_median\_10}', ...
    '\texttt{diff\_mean\_10}', ...
    '\texttt{diff\_median\_10}', ...
    '\texttt{ratio\_mean\_25}', ...
    '\texttt{ratio\_median\_25}', ...
    '\texttt{diff\_mean\_25}', ...
    '\texttt{diff\_median\_25}', ...
    ... ela_distribution
    '$\overline{\phantom{AAAA}}$\qquad\quad\ \texttt{skewness}', ...
    '$\Phi_\textrm{y-D}$ \qquad\quad\texttt{kurtosis}', ...
    '\texttt{num\_of\_peaks}', ...
    ... ela_levelset
    '$\overline{\phantom{AAAA}}$\qquad\texttt{mmce\_lda\_10}', ...
    '\texttt{mmce\_qda\_10}', ...
    '\texttt{mmce\_mda\_10}', ...
    '\texttt{lda\_qda\_10}', ...
    '\texttt{lda\_mda\_10}', ...
    '\texttt{qda\_mda\_10}', ...
    '\texttt{mmce\_lda\_25}', ...
    '\texttt{mmce\_qda\_25}', ...
    '$\Phi_\textrm{Lvl}$ \qquad\texttt{mmce\_mda\_25}', ...
    '\texttt{lda\_qda\_25}', ...
    '\texttt{lda\_mda\_25}', ...
    '\texttt{qda\_mda\_25}', ...
    '\texttt{mmce\_lda\_50}', ...
    '\texttt{mmce\_qda\_50}', ...
    '\texttt{mmce\_mda\_50}', ...
    '\texttt{lda\_qda\_50}', ...
    '\texttt{lda\_mda\_50}', ...
    '\texttt{qda\_mda\_50}', ...
    ... ela_metamodel
    '$\overline{\phantom{AAAA}}$\qquad\ \ \texttt{l\_s\_adj\_r2}', ... 'lin_simple_adj_r2', ...
    '\texttt{l\_s\_c\_min}', ... 'lin_simple_coef_min', ...
    '\texttt{l\_s\_c\_max}', ... 'lin_simple_coef_max', ...
    '$\Phi_\textrm{MM}$\quad\,\texttt{l\_s\_c\_max\_min}', ... 'lin_simple_coef_max_by_min', ...
    '\texttt{l\_w\_i\_adj\_r2}', ... 'lin_w_interact_adj_r2', ...
    '\texttt{q\_s\_adj\_r2}', ... 'quad_simple_adj_r2', ...
    '\texttt{q\_s\_cond}', ... 'quad_simple_cond', ...
    '\texttt{q\_w\_i\_adj\_r2}', ... 'quad_w_interact_adj_r2', ...
    ... infocontent
    '$\overline{\phantom{AAAA}}$\qquad\qquad\quad\texttt{h\_max}', ...
    '\texttt{eps\_s}', ...
    '$\Phi_\textrm{Inf}$\qquad\qquad\ \,\texttt{eps\_max}', ...
    '\texttt{m0}', ...
    '\texttt{eps\_ratio}', ...
    ... nearest_better
    '$\overline{\phantom{AAAA}}$\quad\ \,\texttt{nb\_std\_ratio}', ...
    '\texttt{nb\_mean\_ratio}', ...
    '$\Phi_\textrm{NBC}$ \qquad\quad\phantom{,} \texttt{nb\_cor}', ...
    '\texttt{dist\_ratio}', ...
    '\texttt{nb\_fitness\_cor}', ...
  };
% create same lengths
% imgFeatNames = cellfun(@(x) sprintf('%30s', x), imgFeatNames, 'UniformOutput', false);

% draw coeffs without colorbar
han{1} = figure('Units', 'centimeters', ...
                'Position', [1 1 sizeX sizeY], ...
                'PaperSize', [sizeX + 2, sizeY + 2]);

% adjust data
plotData = ks_sign_2(featNonId, :)';
% non-significant data
inSignifId = plotData > 0.05;
% data closest to significance level
secondSignifVal = 1e-5;
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
ax.XTick = 1:numel(featNonId);
ax.XTickLabel = imgFeatNames;
ax.YTick = 1:(3*nModel);
ax.YTickLabel = modelSetNames;
ax.XTickLabelRotation = labelRot;
ax.TickLabelInterpreter = 'latex';
grid off
% set(ax, 'visible', 'off')
% set(findall(ax, 'type', 'text'), 'visible', 'on')

% change insignificant data to red
han{1}.Colormap(end, :) = [1 0 0];

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
