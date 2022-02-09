%% Surrogate model selection plots
% Script for analysing GP covariance function dependence on data sampled 
% from DTS-CMA-ES run over the noiseless part of the BBOB framework.
% Script also creates graphs and tables.
%

%% 

% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'smsp_data.mat');
if isfile(tmpFName)
  load(tmpFName);
else


%%

% folders for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ecml2021journal'); %'smsppaper');
supplementaryFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ecml2021journal_supp');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
plotSuppResultsFolder = fullfile(supplementaryFolder, 'images');
tableSuppFolder = fullfile(supplementaryFolder, 'tex');
experimentFolder = fullfile('exp', 'experiments', 'exp_smsp');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);
[~, ~] = mkdir(plotSuppResultsFolder);
[~, ~] = mkdir(tableSuppFolder);
[~, ~] = mkdir(experimentFolder);

% create tables from calculated features and save them

expfolder = fullfile('exp', 'experiments');
exp_id = 'exp_DTSmeta_03';
exp_id_knn = 'exp_DTSmeta_03_knn';

dims = [2, 3, 5, 10, 20];
tssList = {'full', 'nearest', 'knn'};
nTSS = numel(tssList);
fullId = find(strcmp(tssList, 'full'));
nearId = find(strcmp(tssList, 'nearest'));
knnId  = find(strcmp(tssList, 'knn'));

exp_output = cellfun(@(x) fullfile(experimentFolder, x), tssList, 'Uni', false);
exp_reduced_output = cellfun(@(x) fullfile(experimentFolder, x, 'sampled_mfts_reduced'), tssList, 'Uni', false);
% make directories
for tss = 1:nTSS
  [~, ~] = mkdir(exp_output{tss});
  [~, ~] = mkdir(exp_reduced_output{tss});
end

exp_meta_minmax   = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_minmax.mat'),   exp_output, 'Uni', false);
exp_meta_inf      = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_inf.mat'),      exp_output, 'Uni', false);
exp_meta_output   = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_res.mat'),      exp_output, 'Uni', false);
exp_meta_quantile = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_quantile.mat'), exp_output, 'Uni', false);
exp_meta_stats    = cellfun(@(x) fullfile(x, 'exp_DTSmeta_03_stats.mat'),    exp_output, 'Uni', false);
exp_meta_dimgen   = fullfile(experimentFolder, 'exp_DTSmeta_03_dimgen.mat');
printScriptMess = false;

% resulting files
exp_smsp_dimension_test = cellfun(@(x) fullfile(x, 'exp_smsp_dimension_test.mat'), exp_output, 'Uni', false);
exp_smsp_corr_test      = cellfun(@(x) fullfile(x, 'exp_smsp_corr_test.mat'),      exp_output, 'Uni', false);
exp_smsp_corr_cluster   = cellfun(@(x) fullfile(x, 'exp_smsp_corr_cluster.mat'),   exp_output, 'Uni', false);
exp_smsp_nan_anal       = cellfun(@(x) fullfile(x, 'exp_smsp_nan_anal.mat'),       exp_output, 'Uni', false);
exp_smsp_quant_anal     = cellfun(@(x) fullfile(x, 'exp_smsp_quant_anal.mat'),     exp_output, 'Uni', false);

exp_smsp_feat_table     = cellfun(@(x) fullfile(tableFolder, ['feat_', x]), tssList, 'Uni', false);
exp_smsp_corr_dim_table = cellfun(@(x) fullfile(tableFolder, ['corr_dim_', x]), tssList, 'Uni', false);

% file lists
fileList = searchFile(fullfile(expfolder, exp_id, 'dataset', 'sampled_metafeatures'), '*.mat*');
nResFiles = numel(fileList);
fileList_knn = searchFile(fullfile(expfolder, exp_id_knn, 'dataset', 'sampled_metafeatures'), '*.mat*');
nResFiles_knn = numel(fileList_knn);

% set random generator seed for replicability
rng = 42;

%%

% Metafeature names

% run experiment definition file to gain settings (ensure first that
% modelTestSets command is disabled)
eval(exp_id)

mfts_group.basic = {'dim', 'observations', 'lower_min', 'lower_max', ...
              'upper_min', 'upper_max', 'objective_min', 'objective_max', ...
              'blocks_min', 'blocks_max', 'cells_total', 'cells_filled', ...
              'minimize_fun'};
mfts_group.cmaes = {'cma_generation', 'cma_step_size', 'cma_restart', ...
              'cma_mean_dist', 'cma_evopath_c_norm', 'cma_evopath_s_norm', ...
              'cma_lik'};
mfts_group.dispersion = {'ratio_mean_02', 'ratio_median_02', ...
                         'diff_mean_02', 'diff_median_02', ...
                         'ratio_mean_05', 'ratio_median_05', ...
                         'diff_mean_05', 'diff_median_05', ...
                         'ratio_mean_10', 'ratio_median_10', ...
                         'diff_mean_10', 'diff_median_10', ...
                         'ratio_mean_25', 'ratio_median_25', ...
                         'diff_mean_25', 'diff_median_25'};
mfts_group.ela_distribution = {'skewness', 'kurtosis', 'number_of_peaks'};
mfts_group.ela_levelset = {'mmce_lda_10', 'mmce_qda_10', 'mmce_mda_10', ...
                           'lda_qda_10', 'lda_mda_10', 'qda_mda_10', ...
                           'mmce_lda_25', 'mmce_qda_25', 'mmce_mda_25', ...
                           'lda_qda_25', 'lda_mda_25', 'qda_mda_25', ...
                           'mmce_lda_50', 'mmce_qda_50', 'mmce_mda_50', ...
                           'lda_qda_50', 'lda_mda_50', 'qda_mda_50'};
mfts_group.ela_metamodel = {'lin_simple_adj_r2', 'lin_simple_intercept', ...
                           'lin_simple_coef_min', 'lin_simple_coef_max', ...
                           'lin_simple_coef_max_by_min', ...
                           'lin_w_interact_adj_r2', 'quad_simple_adj_r2', ...
                           'quad_simple_cond', 'quad_w_interact_adj_r2'};
mfts_group.infocontent = {'h_max', 'eps_s', 'eps_max', 'm0', 'eps_ratio'};
mfts_group.nearest_better = {'nb_std_ratio', 'nb_mean_ratio', 'nb_cor', ...
                             'dist_ratio', 'nb_fitness_cor'};

% mfts name generation cycle
mfts_names = {};
for mi = 1:numel(opts.mfts_settings.MetaInput)
  for uf = opts.mfts_settings.useFeat{mi}
    for mn = 1:numel(mfts_group.(opts.mfts_settings.features{uf}))
      mfts_names{end+1, 1} = strjoin({opts.mfts_settings.MetaInput{mi}, ...
                                      opts.mfts_settings.TransData{mi}, ...
                                      opts.mfts_settings.features{uf}, ...
                                      mfts_group.(opts.mfts_settings.features{uf}){mn}}, ...
                                     '_');
    end
  end
end

% identify last archive feature
lastArchFtId = find(cellfun(@(x) contains(x, 'archive'), mfts_names), 1, 'last');

% clear unnecessary variables
clear mi mn modelOptions uf

%% 

% Removing metafeatures
% Remove user-defined, constant, NaN, and linearly dependent metafeatures.

% Process metafeatures

if printScriptMess
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Processing metafeatures\n'], fix(clock))
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
% 14: ela_metamodel.lin_simple_intercept

rem_mfts = {};
for mi = 1:numel(opts.mfts_settings.MetaInput)
  mss = strjoin({opts.mfts_settings.MetaInput{mi}, ...
                 opts.mfts_settings.TransData{mi}, ''}, '_');
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
% also remove metafeatures identical for all mfts sets except the first
% group
rem2_mfts = {};
for mi = 2:numel(opts.mfts_settings.MetaInput)
  mss = strjoin({opts.mfts_settings.MetaInput{mi}, ...
                 opts.mfts_settings.TransData{mi}, ''}, '_');
  rem2_mfts = [rem2_mfts, {...
      [mss, 'basic_dim'], ...
      [mss, 'cmaes_cma_evopath_c_norm'], ...
      [mss, 'cmaes_cma_evopath_s_norm'], ...
      [mss, 'cmaes_cma_generation'], ...
      [mss, 'cmaes_cma_restart'], ...
      [mss, 'cmaes_cma_step_size'] ...
    }];
end
% remove infocontent metafeatures dependent on set M for test sets -> M is
% not computed correctly if y-values are missing
rem3_mfts = {};
for mi = 1:numel(opts.mfts_settings.MetaInput)
  if contains(opts.mfts_settings.MetaInput{mi}, 'test')
    mss = strjoin({opts.mfts_settings.MetaInput{mi}, ...
                   opts.mfts_settings.TransData{mi}, ''}, '_');
    rem3_mfts = [rem3_mfts, {...
        [mss, 'infocontent_m0'], ...
        [mss, 'infocontent_eps_ratio'], ...
      }];
  end
end
rem_mfts = [rem_mfts, rem2_mfts, rem3_mfts];

fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
         'Removing user-defined metafeatures:\n'], fix(clock))
removeId = false(numel(mfts_names), 1);
for m = 1:numel(rem_mfts)
  mId = strcmp(mfts_names, rem_mfts{m});
  fprintf('%s\n', rem_mfts{m})
  removeId(mId) = true;
end
% mfts = mfts(:, ~mId);
mfts_names = mfts_names(~removeId);
mftsIds = ~removeId;
% vars = vars(~removeId, :);
% means = means(~removeId, :);

% rename archive identical columns (see above)
identicalLabel = strjoin({opts.mfts_settings.MetaInput{1}, ...
                 opts.mfts_settings.TransData{1}, ''}, '_');
mfts_names{strcmp(mfts_names, [identicalLabel, 'basic_dim'])} = 'dim';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_evopath_c_norm'])} = 'cmaes_evopath_c_norm';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_evopath_s_norm'])} = 'cmaes_evopath_s_norm';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_generation'])}     = 'cmaes_generation';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_restart'])}        = 'cmaes_restart';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_step_size'])}      = 'cmaes_step_size';
% rename observations columns
for m = 1:numel(mfts_names)
  if contains(mfts_names{m}, 'observations')
    mfts_name_parts = strsplit(mfts_names{m}, '_');
    mfts_names{m} = strjoin({mfts_name_parts{1}, 'obs'}, '_');
  end
end

% create mfts table, where only user defined metafeatures are not present
% full_mfts = mfts;

% metafeature names for printing
mfts_names_prt = strrep(mfts_names, '_', '\_');

% identify last archive feature after reduction
lastArchFtId_red = find(cellfun(@(x) contains(x, 'archive'), mfts_names), 1, 'last');

clear rem_mfts rem2_mfts rem3_mfts m ms mId mss

%% 

% reduce results to feature values, names, and basic identifiers

% list of ids to reduce (exclude n.7)
redIds = [1:6, 8:9];
errListFile = '/tmp/errFileList';
errListFile_knn = '/tmp/errFileList_knn';

% reduce TSS nearest
for fl = 1:nResFiles
  [~, actFile] = fileparts(fileList{fl});
  reducedFilename_full = fullfile(exp_reduced_output{fullId}, [actFile, '_red.mat']);
  reducedFilename_near = fullfile(exp_reduced_output{nearId}, [actFile, '_red.mat']);
  id_str = regexp(actFile, '_id(\d)_', 'tokens');
  id = str2double(id_str{1});
  % create reduced TSS nearest and TSS full results
  if (~isfile(reducedFilename_full) || ~isfile(reducedFilename_near)) ...
      && ismember(id, redIds)
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Saving %s\n'], fix(clock), reducedFilename_full)
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Saving %s\n'], fix(clock), reducedFilename_near)
    end
    try
      S = load(fileList{fl});
    catch err
      warning(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'File %s loading failed with the following error: %s'], fix(clock), fileList{fl}, err.message)
      FID = fopen(errListFile, 'a');
      fprintf(FID, '%s\n', fileList{fl});
      fclose(FID);
      clear FID
      continue
    end
    dim = S.dim;
    fun = S.fun;
    inst = S.inst;
    % TSS full
    mfts_values = S.values(~removeId(1:lastArchFtId), :, :);
    mfts_names_red = mfts_names(1:lastArchFtId_red);
    save(reducedFilename_full, 'dim', 'fun', 'inst', 'id', 'mfts_names_red', 'mfts_values')
    % TSS nearest
    mfts_values = S.values(~removeId(lastArchFtId+1:end), :, :);
    mfts_names_red = mfts_names(lastArchFtId_red+1:end);
    save(reducedFilename_near, 'dim', 'fun', 'inst', 'id', 'mfts_names_red', 'mfts_values')
  elseif isfile(reducedFilename_full)
%     S = load(fileList{fl});
  end
end

for fl = 1:nResFiles_knn
  [~, actFile] = fileparts(fileList_knn{fl});
  % create reduced TSS knn results
  fileName_knn = fullfile(expfolder, exp_id_knn, 'dataset', 'sampled_metafeatures', [actFile, '.mat']);
  reducedFilename_knn = fullfile(exp_reduced_output{knnId}, [actFile, '_red.mat']);
  id_str = regexp(actFile, '_id(\d)_', 'tokens');
  id = str2double(id_str{1});
  if ~isfile(reducedFilename_knn) && ismember(id, redIds)
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Saving %s\n'], fix(clock), reducedFilename_knn)
    end
    try
      S_knn = load(fileList_knn{fl});
    catch err
      warning(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'File %s loading failed with the following error: %s'], fix(clock), fileList_knn{fl}, err.message)
      FID = fopen(errListFile_knn, 'a');
      fprintf(FID, '%s\n', fileList_knn{fl});
      fclose(FID);
      clear FID
      continue
    end
    dim = S_knn.dim;
    fun = S_knn.fun;
    inst = S_knn.inst;
    % TSS knn
    mfts_values = S_knn.values(~removeId(lastArchFtId+1:end), :, :);
    mfts_names_red = mfts_names(lastArchFtId_red+1:end);
    save(reducedFilename_knn, 'dim', 'fun', 'inst', 'id', 'mfts_names_red', 'mfts_values')
  end
end

clear actFile filename_knn fl dim fun inst id id_str mfts_values S S_knn ...
      reducedFilename_full reducedFilename_nearest reducedFilename_knn

%% 

% feature distribution - is it necessary - code is not effective for large
% number of files

q_levels = 0:0.01:1;

% % TSS loop
% for tss = 1:nTSS
%   if ~isfile(exp_meta_quantile{tss})
%     % create list of files with reduced data
%     redFileList = searchFile(exp_reduced_output{tss}, '*.mat*');
%     % load first file to get metafeature names
%     S = load(redFileList{1}, 'mfts_names_red');
%     mfts_names_red = S.mfts_names_red;
%     for ft = 1:numel(mfts_names_red)
%       fprintf('Distribution TSS %s (%3d/%3d): %s\n', tssList{tss}, ...
%               ft, numel(mfts_names_red), mfts_names_red{ft})
%
%       actual_ft = [];
%       for fl = 1 : numel(redFileList)
%         S = load(redFileList{fl}, 'mfts_values');
%         actual_ft = [actual_ft, shiftdim(S.mfts_values(ft, :, :), 1)];
%       end
%       q_ft(ft, :) = quantile(actual_ft(:), q_levels);
%     end
%     save(exp_meta_quantile{tss}, 'q_ft', 'mfts_names_red')
%   end
% end

clear actual_ft fl ft q_levels S

%% 

% feature mins and maxs

% TSS loop
for tss = 1:nTSS
  % create list of files with reduced data
  redFileList = searchFile(exp_reduced_output{tss}, '*.mat*');

  if isfile(exp_meta_minmax{tss})
%     MM = load(exp_meta_minmax{tss});
%     ft_max = MM.ft_max;
%     ft_min = MM.ft_min;
%     ft_min_ninf = MM.ft_min_ninf;
%     ft_max_ninf = MM.ft_max_ninf;
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'File with minimums and maximums for TSS %s already exists.\n'], ...
              fix(clock), tssList{tss})
  else
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Calculating feature minimums and maximums for TSS %s\n'], ...
              fix(clock), tssList{tss})
    end
    S = load(redFileList{1}, 'mfts_names_red');
    mfts_names_red = S.mfts_names_red;
    nMftsAct = numel(mfts_names_red);

    % init extreme values
    ft_min = NaN(nMftsAct, 1);
    ft_max = NaN(nMftsAct, 1);
    ft_min_ninf = ft_min;
    ft_max_ninf = ft_max;
    S_ft_min_ninf = NaN(size(ft_min_ninf));
    S_ft_max_ninf = NaN(size(ft_max_ninf));

    % compare with the rest of results
    for fl = 1:numel(redFileList)
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'MinMax TSS %s (%4d/%4d): %s\n'], ...
              fix(clock), tssList{tss}, ...
              fl, numel(redFileList), redFileList{fl})
      S = load(redFileList{fl}, 'mfts_values');
      ft_min = min(ft_min, min(min(S.mfts_values, [], 3), [], 2));
      ft_max = max(ft_max, max(max(S.mfts_values, [], 3), [], 2));
      % find non-inf min values
      for ft = find(isinf(ft_min))
        actMfts = S.mfts_values(ft, :, :);
        act_ft_min_ninf = min(actMfts(~isinf(actMfts)));
        if ~isempty(act_ft_min_ninf)
          S_ft_min_ninf(ft, 1) = act_ft_min_ninf;
        end
      end
      % find non-inf max values
      for ft = find(isinf(ft_max))
        actMfts = S.mfts_values(ft, :, :);
        act_ft_max_ninf = max(actMfts(~isinf(actMfts)));
        if ~isempty(act_ft_max_ninf)
          S_ft_max_ninf(ft, 1) = act_ft_max_ninf;
        end
      end
      % compare to overall 2nd min and max
      ft_min_ninf = min(ft_min_ninf, S_ft_min_ninf);
      ft_max_ninf = max(ft_max_ninf, S_ft_max_ninf);
    end

    % replace NaNs in non-inf metafeatures with real min and max
    ft_min_ninf(~isinf(ft_min) & isnan(ft_min_ninf)) = ft_min(~isinf(ft_min));
    ft_max_ninf(~isinf(ft_max) & isnan(ft_max_ninf)) = ft_max(~isinf(ft_max));

    % save inf and non-inf min and max
    save(exp_meta_minmax{tss}, 'ft_min', 'ft_max', 'ft_min_ninf', 'ft_max_ninf', ...
                               'mfts_names_red')
  end
end

clear actMfts act_ft_max act_ft_min fl ft ft_min ft_max ft_min_ninf ...
      ft_max_ninf mfts_names_red MM redFileList S S_ft_max_ninf ...
      S_ft_min_ninf tss

%%

% feature infs

% TSS loop
for tss = 1:nTSS
  % create list of files with reduced data
  redFileList = searchFile(exp_reduced_output{tss}, '*.mat*');

  if isfile(exp_meta_inf{tss})
%     MM = load(exp_meta_inf);
%     inf_plus = MM.inf_plus;
%     inf_min = MM.inf_min;
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'File with infs for TSS %s already exists.\n'], ...
              fix(clock), tssList{tss})
  else
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Calculating feature infs\n'], fix(clock))
    end
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
             'Inf TSS %s (%3d/%3d): %s\n'], fix(clock), tssList{tss}, ...
            1, numel(redFileList), redFileList{1})

    S = load(redFileList{1}, 'mfts_values', 'mfts_names_red');
    mfts_names_red = S.mfts_names_red;
    % calculate number of infs
    inf_plus = sum(sum(isinf(S.mfts_values) & S.mfts_values > 0, 3), 2);
    inf_min = sum(sum(isinf(S.mfts_values) & S.mfts_values < 0, 3), 2);

    % compare with the rest of results
    for fl = 2:numel(redFileList)
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Inf TSS %s (%3d/%3d): %s\n'], fix(clock), tssList{tss}, ...
              fl, numel(redFileList), redFileList{fl})
      S = load(redFileList{fl}, 'mfts_values');
      % calculate number of infs
      inf_plus = inf_plus + sum(sum(isinf(S.mfts_values) & S.mfts_values > 0, 3), 2);
      inf_min = inf_min + sum(sum(isinf(S.mfts_values) & S.mfts_values < 0, 3), 2);
    end

    % save min and max
    save(exp_meta_inf{tss}, 'inf_plus', 'inf_min', 'mfts_names_red')
  end
end

clear fl inf_plus inf_min MM redFileList S tss

%%

% feature quantiles

quantile_vals = [0.01, 0.99];

% TSS loop
for tss = 1:nTSS
  % create list of files with reduced data
  redFileList = searchFile(exp_reduced_output{tss}, '*.mat*');

  if isfile(exp_meta_quantile{tss})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'File with quantiles for TSS %s already exists.\n'], ...
              fix(clock), tssList{tss})
  else
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Calculating feature quantiles\n'], fix(clock))
    end

    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
             'Quantiles TSS %s (%3d/%3d): %s\n'], fix(clock), tssList{tss}, ...
            1, numel(redFileList), redFileList{1})

    S = load(redFileList{1}, 'mfts_values', 'mfts_names_red');
    mfts_names_red = S.mfts_names_red;

    nfts = numel(mfts_names_red);
    kft = 25; % k features per iteration

    ft_vals = NaN(kft, size(S.mfts_values, 2), size(S.mfts_values, 3), numel(redFileList));
    quantiles = NaN(numel(mfts_names_red), 2);

    % feature group cycle (groups according to kft value)
    for ftg = 1:ceil(nfts/kft)
      for fl = 1:numel(redFileList)
        fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Loading values for quantiles TSS %s feature group %d/%d (%3d/%3d): %s\n'], fix(clock), tssList{tss}, ...
              ftg, ceil(nfts/kft), fl, numel(redFileList), redFileList{fl})
        S = load(redFileList{fl}, 'mfts_values');
        % check final feature index
        lastFt = min(ftg*kft, nfts);
        lastFtg = mod(lastFt, kft);
        if lastFtg == 0
          lastFtg = kft;
        end
        ft_vals(1:lastFtg, :, :, fl) = S.mfts_values(((ftg-1)*kft+1):lastFt, :, :);
      end
      for ft = 1:lastFtg
        fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Quantiles calculation TSS %s feature group %d/%d feature %d/%d\n'], fix(clock), tssList{tss}, ...
              ftg, ceil(nfts/kft), ft, lastFtg)
        act_ft_vals = ft_vals(ft, :, :, :);
        % remove Infs
        act_ft_vals = act_ft_vals(~isinf(act_ft_vals));
        quantiles((ftg-1)*kft+ft, :) = quantile(act_ft_vals(:), quantile_vals);
      end
      % save quantiles
      save(exp_meta_quantile{tss}, 'quantiles', 'quantile_vals', 'mfts_names_red')
    end
  end
end

clear act_ft_vals fl ft ft_vals mfts_names_red quantiles quantile_vals redFileList S tss

%% load and normalize results

% quantile settings
% reliability of all 42 future-clustered feature representatives
demanded_reliability = [0.9, 0.99];
nClusters = 42;
quantile_to_calc = @(reliability) [  (1-reliability.^(1/nClusters))/2, ...
                                   1-(1-reliability.^(1/nClusters))/2];

% TSS loop
for tss = 1:nTSS
  % create list of files with reduced data
  redFileList = searchFile(exp_reduced_output{tss}, '*.mat*');

  if isfile(exp_meta_minmax{tss})
    load(exp_meta_minmax{tss}, 'ft_min', 'ft_max', 'ft_min_ninf', 'ft_max_ninf', 'mfts_names_red');
  else
    error('There is %s missing.', exp_meta_minmax{tss})
  end

  if isfile(exp_meta_quantile{tss})
    S_q = load(exp_meta_quantile{tss}, 'quantiles', 'quantile_vals');
  else
    error('There is %s missing.', exp_meta_quantile{tss})
  end

  if isfile(exp_meta_output{tss})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'File with normalized results for TSS %s already exists.\n'], ...
          fix(clock), tssList{tss})
  else
    if printScriptMess
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Loading metalearing results\n'], fix(clock))
    end

    % load first file to find sizes
    S = load(redFileList{1}, 'mfts_values');
    nSamples = size(S.mfts_values, 3);
    nFiles = numel(redFileList);
    nMfts = size(S.mfts_values, 1);

    % init
    observations = NaN(2, nFiles*nSamples); % expected range [0, 5000]
    if tss == fullId
      generations = NaN(1, nFiles*nSamples); % expected range [0, 5000]
      dimensions  = NaN(1, nFiles*nSamples); % expected range [2, 20]
    end
    vars     = NaN(nMfts, nFiles*nSamples);
    means    = NaN(nMfts, nFiles*nSamples);
    medians  = NaN(nMfts, nFiles*nSamples);
    nans     = NaN(nMfts, nFiles*nSamples); % expected range [0, 100]
    inf_plus = NaN(nMfts, nFiles*nSamples); % expected range [0, 100]
    inf_mins = NaN(nMfts, nFiles*nSamples); % expected range [0, 100]
    norm_val = NaN(size(S.mfts_values));
    
    quantiles_90 = NaN(nMfts, nFiles*nSamples, 2);
    quantiles_99 = NaN(nMfts, nFiles*nSamples, 2);

    for fl = 1:nFiles
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Normalization TSS %s (%4d/%4d): %s\n'], ...
              fix(clock), tssList{tss}, ...
              fl, nFiles, redFileList{fl})
      S = load(redFileList{fl}, 'mfts_values');

      % normalize to [0,1]
      for ft = 1:size(S.mfts_values, 1)
        act_mfts_values = S.mfts_values(ft, :, :);

        if (ft_max_ninf(ft)-ft_min_ninf(ft)) == 0
          % constant feature
          norm_val(ft, :, :) = 0.5*ones(size(act_mfts_values));
        else
          % linear transformation to [0, 1]
          % norm_val(ft, :, :) = (act_mfts_values - ft_min_ninf(ft))/(ft_max_ninf(ft)-ft_min_ninf(ft));
          % sigmoid transformation to [0, 1] tuned according to S_q.quantiles
          if S_q.quantile_vals(1) == 1-S_q.quantile_vals(2)
            % symmetric quantiles
            norm_val(ft, :, :) = 1./(1+exp(log((1-S_q.quantile_vals(1))/S_q.quantile_vals(1))* ...
              (2*act_mfts_values - S_q.quantiles(ft, 1)-S_q.quantiles(ft, 2)) / ...
              (S_q.quantiles(ft, 1)-S_q.quantiles(ft, 2))));
          else
            % arbitrary quantiles
            n1 = S_q.quantile_vals(1);
            n2 = S_q.quantile_vals(2);
            k = log(n2*(1-n1)/(n1*(1-n2)))/(S_q.quantiles(ft, 2)-S_q.quantiles(ft, 1));
            x0 = (log((1-n1)*(1-n2)/(n1*n2))*(S_q.quantiles(ft, 2)-S_q.quantiles(ft, 1))/ ...
                  log(n2*(1-n1)/(n1*(1-n2)))+(S_q.quantiles(ft, 2)+S_q.quantiles(ft, 1)))/2;
            norm_val(ft, :, :) = 1./(1+exp(-k*(act_mfts_values-x0)));
          end
        end
      end

      % row id for collection
      colId = (fl-1)*nSamples+1 : fl*nSamples;
      % collect generations, numbers of points, and dimensions
      observations(:, colId) = S.mfts_values(contains(mfts_names_red, 'obs'), :, 1);
      if tss == fullId
        generations(:, colId) = S.mfts_values(contains(mfts_names_red, 'generation'),   :, 1);
        dimensions(:, colId)  = S.mfts_values(contains(mfts_names_red, 'dim'),    :, 1);
      end

      % collect nans and infs
      nans(:, colId) = sum(isnan(norm_val), 3);
      inf_plus(:, colId) = sum(isinf(norm_val) & norm_val > 0, 3);
      inf_mins(:, colId) = sum(isinf(norm_val) & norm_val < 0, 3);

      % calculate statistics from non-inf and non-nan values
      norm_val(isinf(norm_val)) = NaN;
      vars(:, colId)  = nanvar(norm_val, 0, 3);
      means(:, colId) = nanmean(norm_val, 3);
      medians(:, colId) = nanmedian(norm_val, 3);
      quantiles_90(:, colId, :) = quantile(norm_val, quantile_to_calc(demanded_reliability(1)), 3);
      quantiles_99(:, colId, :) = quantile(norm_val, quantile_to_calc(demanded_reliability(2)), 3);

      % save results after each 100 files
      if mod(fl, 100) == 0
        save(exp_meta_output{tss}, 'observations', 'nans', ...
          'inf_plus', 'inf_mins', 'vars', 'means', 'medians', ...
          'quantiles_90', 'quantiles_99', ...
          'mfts_names_red', 'quantile_to_calc')
        % save dimensions and generations
        if tss == fullId
          save(exp_meta_dimgen, 'dimensions', 'generations')
        end
      end
    end

    % save dimensions and generations
    if tss == fullId
      save(exp_meta_dimgen, 'dimensions', 'generations')
    end
    S = load(exp_meta_dimgen, 'dimensions', 'generations');
    dimensions = S.dimensions;
    generations = S.generations;
    % save overall results
    save(exp_meta_output{tss}, 'observations', 'generations', 'dimensions', 'nans', ...
      'inf_plus', 'inf_mins', 'vars', 'means', 'medians', ...
      'quantiles_90', 'quantiles_99', ...
      'mfts_names_red', 'quantile_to_calc')
  %   save(exp_meta_output, 'nans', 'inf_plus', 'inf_mins', 'vars', 'means', 'medians')
  end
end

clear actualFolder articleFolder fl ft i infBound logscale norm_val S act_ft_min act_ft_max act_mfts_values rowId

%%

% calculate simple stats
for tss = 1:nTSS
  
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Simple stats of TSS %s\n'], ...
              fix(clock), tssList{tss})

  % load results separately due to memory requirements
  if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'nans', 'inf_plus', 'inf_mins', ...
                                   'mfts_names_red');
  else
    error('There is %s missing.', exp_meta_output{tss})
  end
  mfts_names_red = S.mfts_names_red;
  nCases = size(S.nans, 2);

  nan_perc = sum(S.nans, 2) / nCases;
  inf_perc = sum(S.inf_mins+S.inf_plus, 2) / nCases;
  
  % clear nan and inf results and load min max results
  clear S
  S = load(exp_meta_output{tss}, 'quantiles_99');

  qDiff = S.quantiles_99(:, :, 2) - S.quantiles_99(:, :, 1);
  maxmin_perc = (sum(~(qDiff > 0.05), 2) - sum(isnan(qDiff), 2)) ./ sum(~isnan(qDiff), 2);

  % mark features with high NaN percentage
  mftsLowNaN = nan_perc < 25;
  mftsLowErr = maxmin_perc > 0.9;

  maxmin_perc_rel = maxmin_perc(mftsLowNaN & mftsLowErr);
  % save stats
  save(exp_meta_stats{tss}, 'inf_perc', 'nan_perc', 'maxmin_perc', ...
                            'mftsLowNaN', 'mftsLowErr', ...
                            'maxmin_perc_rel', 'mfts_names_red')
end

clear('inf_perc', 'nan_perc', 'maxmin_perc', 'maxmin_perc_rel', ...
      'mftsLowNaN', 'mftsLowErr', 'nCases', 'qDiff',  'S', 'tss')

%%

% prepare vars and means for observation and density point of view

% if isfile(exp_meta_output)
%   load(exp_meta_output)
% else
%   error('File %s is missing!', exp_meta_output)
% end
% 
% % overall results
% if printScriptMess
%   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
%            'Overall results\n'], fix(clock))
% end
% overall_var  = nanmedian(vars, 2);
% overall_mean = nanmedian(means, 2);
% 
% % results according to dimension
% if printScriptMess
%   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
%            'Dimension results\n'], fix(clock))
% end
% 
% for d = 1:numel(dims)
%   dim_var(:, d)  = nanmedian(vars(:, dims(d) == dimensions), 2);
%   dim_mean(:, d) = nanmedian(means(:, dims(d) == dimensions), 2);
% end
% 
% % number of xaxis quantiles
% nQuant = 1000;
% nMfts = numel(mfts_names);
% % plot data quantiles
% disp_quant = [0.05, 0.5, 0.95];
% nDispQuant = numel(disp_quant);
% 
% % init variables
% obs_vars = NaN(nMfts, nQuant, nDispQuant);
% obs_means = NaN(nMfts, nQuant, nDispQuant);
% obs_nans = NaN(nMfts, nQuant, nDispQuant);
% obs_inf_plus = NaN(nMfts, nQuant, nDispQuant);
% obs_inf_mins = NaN(nMfts, nQuant, nDispQuant);
% obs_dims = false(size(observations, 1), nQuant, numel(dims));
% 
% dens_vars = NaN(nMfts, nQuant, nDispQuant);
% dens_means = NaN(nMfts, nQuant, nDispQuant);
% dens_nans = NaN(nMfts, nQuant, nDispQuant);
% dens_inf_plus = NaN(nMfts, nQuant, nDispQuant);
% dens_inf_mins = NaN(nMfts, nQuant, nDispQuant);
% dens_dims = false(size(observations, 1), nQuant, numel(dims));
% 
% if printScriptMess
%   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
%            'Results according to the number of observations\n'], ...
%           fix(clock))
% end
% sampleSets = {'archive_', 'archivetest_', 'train_', 'traintest_'};
% % calculate variance and mean medians according to number of observations 
% % for different sample sets
% noSampleSetId = ~contains(mfts_names, sampleSets);
% nGenFeat = sum(noSampleSetId); % number of sample independent features
% un_observations = cell(1, 3);
% for o = 1:size(observations, 1)
%   un_observations{o} = unique(observations(o, :));
%   sampleSetId = contains(mfts_names, sampleSets{o});
%   % calculate densities
%   densities(o, :) = nthroot(observations(o, :), dimensions);
%   % get ids of sorted values
%   [~,  obs_I] = sort(observations(o, :));
%   [~, dens_I] = sort(   densities(o, :));
%   % exclude zero observations
%   obs_I(~isnatural(observations(o, obs_I))) = [];
%   dens_I(~(densities(o, dens_I) > 0)) = [];
%   % create permilles (1000-quantile) data groups
%   obs_gId  = ceil(nQuant*(1:numel(obs_I))  / numel(obs_I));
%   dens_gId = ceil(nQuant*(1:numel(dens_I)) / numel(dens_I));
%   
%   for q = 1:nQuant
%     % observations
%     act_dataId = obs_I(obs_gId == q);
%     obs_vars(sampleSetId, q, :)      = shiftdim(quantile(     vars(sampleSetId, act_dataId), disp_quant, 2), -1);
%     obs_means(sampleSetId, q, :)     = shiftdim(quantile(    means(sampleSetId, act_dataId), disp_quant, 2), -1);
%     obs_nans(sampleSetId, q, :)      = shiftdim(quantile(     nans(sampleSetId, act_dataId), disp_quant, 2), -1);
%     obs_inf_plus(sampleSetId, q, :)  = shiftdim(quantile( inf_plus(sampleSetId, act_dataId), disp_quant, 2), -1);
%     obs_inf_mins(sampleSetId, q, :)  = shiftdim(quantile( inf_mins(sampleSetId, act_dataId), disp_quant, 2), -1);
%     % dimensions present in particular quantile
%     obs_dims(o, q, :) = ismember(dims, dimensions(act_dataId));
%     % observations in particular quantile
%     obs_obs(o, q, :) = median(observations(o, act_dataId));
%     % calculate only means for sample independent features
% %     general_obs_means((1:nGenFeat) + (o-1)*nGenFeat, uo) = ...
% %       median(means(noSampleSetId, un_observations{o}(uo) == observations(o, :)), 2);
% 
%     % densities
%     act_dataId = dens_I(dens_gId == q);
%     dens_vars(sampleSetId, q, :)      = shiftdim(quantile(     vars(sampleSetId, act_dataId), disp_quant, 2), -1);
%     dens_means(sampleSetId, q, :)     = shiftdim(quantile(    means(sampleSetId, act_dataId), disp_quant, 2), -1);
%     dens_nans(sampleSetId, q, :)      = shiftdim(quantile(     nans(sampleSetId, act_dataId), disp_quant, 2), -1);
%     dens_inf_plus(sampleSetId, q, :)  = shiftdim(quantile( inf_plus(sampleSetId, act_dataId), disp_quant, 2), -1);
%     dens_inf_mins(sampleSetId, q, :)  = shiftdim(quantile( inf_mins(sampleSetId, act_dataId), disp_quant, 2), -1);
%     % dimensions present in particular quantile
%     dens_dims(o, q, :) = ismember(dims, dimensions(act_dataId));
%     % densities in particular quantile
%     dens_dens(o, q, :) = median(densities(o, act_dataId));
%   end
% end

clear d o q
% clear large variables saved in exp_meta_output
clear means medians vars nans inf_plus inf_mins

%%

% save temporary variables for further usage
if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Test dependency of feature median on dimension
% Rejecting independence of medians on dimension is equivalent to rejecting
% equality of medians for any pair of dimensions using Bonferroni
% correction on the number of possible pairs.

alpha = 0.05;

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing Friedman test on dimension TSS %s\n'], ...
          fix(clock), tssList{tss})

  if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'dimensions', 'medians', ...
                                   'mfts_names_red');
  else
    error('File %s is missing!', exp_meta_output{tss})
  end
  
  % use all features
  mfts_names_red = S.mfts_names_red;
  medians = S.medians;

  nMfts = numel(mfts_names_red);
  % medtest_meds_p = zeros(nMfts, nchoosek(numel(dims), 2));
  friedman_meds_p = NaN(nMfts, 1);
  stats_meds = cell(nMfts, 1);
  multcomp_meds = cell(nMfts, 1);
  friedman_meds_p_bh = NaN(nMfts, nchoosek(numel(dims), 2));
  wilcoxon_meds_p    = NaN(nMfts, nchoosek(numel(dims), 2));
  wilcoxon_meds_p_bh = NaN(nMfts, nchoosek(numel(dims), 2));
  for m = 1:nMfts
    act_meds = medians(m, :);
  %   j = 0;
  %   for d = 1:numel(dims)
  %     for d2 = d+1 : numel(dims)
  %       % increase counter
  %       j = j+1;
  %       median test on all non-nan median values of a pair of dimensions
  %       medtest_meds_p(m, j) = mediantest(act_meds(~isnan(act_meds) & dims(d) == dimensions), ...
  %                                         act_meds(~isnan(act_meds) & dims(d2) == dimensions));
  %     end
  %   end

  % % display independent features on the level alpha/number of pairs
  % % (Bonferroni)
  %   if all(medtest_meds_p(m, :) > alpha/nchoosek(numel(dims), 2))
  %     fprintf('Not rejecting independence of medians on dimension for %s on (alpha = %0.3f) using Bonferroni correction\n', mfts_names{m}, alpha)
  %   end

    % Friedman's test
    mat_meds = [];
    % find median sizes
    dimSizes = arrayfun(@(x) sum(x == S.dimensions), dims);
    minDimSize = min(dimSizes);
    if any(dimSizes ~= minDimSize)
      warning('The size of results in different dimensions are not the same for %s. Reducing to the lowest size.', tssList{tss})
    end
    % cat medians for individual dimensions (number of data for each
    % dimension must be identical)
    for d = 1:numel(dims)
      mat_meds_act = act_meds(dims(d) == S.dimensions)';
      mat_meds = [mat_meds, mat_meds_act(1:minDimSize)];
    end
    % remove NaN's
    mat_meds(any(isnan(mat_meds), 2), :) = [];
    if size(mat_meds, 1) > 1 && size(mat_meds, 2) > 1
      % Friedman's test on statistical significance of pairwise differences
      % of dimension medians
      [friedman_meds_p(m), ~, stats_meds{m}] = friedman(mat_meds, 1, 'off');
    %   fprintf('%48s Friedman: %f\n', mfts_names{m}, friedman_meds_p(m))
      multcomp_meds{m} = multcompare(stats_meds{m});
      % p-values using Bonferroni-Holm correction on the alpha level
      friedman_meds_p_bh(m, :) = bonfHolm(multcomp_meds{m}(:, end), alpha);
      % Wilcoxon signed rank test
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
               'Performing two-sided Wilcoxon signed rank test on pairs of feature dimensions TSS %s\n'], ...
              fix(clock), tssList{tss})
      j = 0;
      for d = 1:numel(dims)
        for d2 = d+1 : numel(dims)
          % increase counter
          j = j+1;
          meds_pair = [act_meds(dims(d) == S.dimensions); ...
                       act_meds(dims(d2) == S.dimensions)];
          % Wilcoxon test on all non-nan median values of dimension pair
          wilcoxon_meds_p(m, j) = signrank(meds_pair(1, all(~isnan(meds_pair))), ...
                                           meds_pair(2, all(~isnan(meds_pair))));
        end
      end
      % p-values using Bonferroni-Holm correction on the alpha level
      wilcoxon_meds_p_bh(m, :) = bonfHolm(wilcoxon_meds_p(m, :), alpha);
    end

  end

  % save testing results
  save(exp_smsp_dimension_test{tss}, ...
    'alpha', 'dims', 'mfts_names_red', ...
    'friedman_meds_p', 'stats_meds', 'multcomp_meds', 'friedman_meds_p_bh', ...
    'wilcoxon_meds_p', 'wilcoxon_meds_p_bh')
end

% clear large variables saved in exp_meta_output
clear('alpha', 'd', 'm', 'mat_meds', 'meds_pair', ...
  'mfts_names_red', ...
  'friedman_meds_p', 'stats_meds', 'multcomp_meds', 'friedman_meds_p_bh', ...
    'wilcoxon_meds_p', 'wilcoxon_meds_p_bh')

%% Feature dimensional influence clustering
% Split metafeatures to groups according to the rejected pairs of
% dimensions

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing feature dimensional inluence clustering on TSS %s\n'], ...
          fix(clock), tssList{tss})
  if isfile(exp_smsp_dimension_test{tss})
    S = load(exp_smsp_dimension_test{tss}, 'alpha', 'wilcoxon_meds_p_bh');
  else
    error('File %s is missing!', exp_smsp_dimension_test{tss})
  end

  % find unique combinations of rejected hypothesis
%   [friedman_uniq_combs, ~, friedman_uniq_combs_id] = unique(friedman_meds_p_bh > alpha, 'rows');
  [wilcoxon_uniq_combs, ~, wilcoxon_uniq_combs_id] = unique(S.wilcoxon_meds_p_bh > S.alpha, 'rows');

  % add results to already existing ones
  save(exp_smsp_dimension_test{tss}, 'wilcoxon_uniq_combs', 'wilcoxon_uniq_combs_id', '-append')
end

clear('alpha', 'S', 'tss', ...
      'wilcoxon_uniq_combs', 'wilcoxon_uniq_combs_id')

%% Analyse NaN values
% Calculate numbers of NaNs over the number of evaluations for each
% feature. Find the number when the majority of feature values can be
% calculated (99%).

nan_threshold = 0.01;

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing NaN values analysis on TSS %s\n'], ...
          fix(clock), tssList{tss})

  % load results
  if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'observations', ...
             'nans', ...
             'mfts_names_red');
  else
    error('File %s is missing!', exp_meta_output{tss})
  end

  % set keyword according to TSS method
  if tss == fullId
    keyword = 'archive';
  else
    keyword = 'train';
  end

  % find unique observations values
  obsUni{1} = unique(S.observations(1, :));
  obsUni{2} = unique(S.observations(2, :));

  % use all features
  mfts_names_red = S.mfts_names_red;
  
  % init
  nMfts = numel(mfts_names_red);
  nan_break = zeros(nMfts, 1);
  nan_meds = cell(nMfts, 1);
  % feature loop
  for ft = 1:nMfts

    %%
    close
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Processing NaN feature analysis on TSS %s feature %3d/%d\n'], ...
          fix(clock), tssList{tss}, ft, nMfts)
    % find observation row
    if contains(mfts_names_red{ft}, [keyword, 'test'])
      obsId = 2;
    elseif contains(mfts_names_red{ft}, keyword)
      obsId = 1;
    else
      continue
    end

    % gather nan medians for unique observation values
    for i = 1:numel(obsUni{obsId})
      uniObsId = obsUni{obsId}(i) == S.observations(obsId, :);
      nan_meds{ft}(i) = median(S.nans(ft, uniObsId));
    end

    % find the breaking point when the number of point is sufficient for
    % 99% of data
    nan_breakId = 1;
    nan_break(ft) = obsUni{obsId}(nan_breakId);
    while sum(S.nans(ft, S.observations(obsId, :) >= nan_break(ft)))/size(S.observations(obsId, :), 2)/100 > nan_threshold
      nan_breakId = nan_breakId + 1;
      nan_break(ft) = obsUni{obsId}(nan_breakId);
    end

    % plot resulting data
    figure()
    plot(obsUni{obsId}, nan_meds{ft})
    hold on
    scatter(nan_break(ft), 0, 200, 'r', '+')
    title(['TSS ', num2str(tss), ': ', strrep(mfts_names_red{ft}, '_', '\_')])
    legend('% of NaN', ['threshold: ', num2str(nan_break(ft)), ' points'])
    hold off
  end

  % save results
  save(exp_smsp_nan_anal{tss}, 'nan_break', 'nan_meds', 'nan_threshold', ...
                               'mfts_names_red', 'obsUni')
end

clear('ft', 'i', 'mfts_names_red', 'nan_break', 'nan_breakId', ...
      'nan_meds', 'nMfts', 'obsId', 'obsUni', 'tss')

close

%% Reliability threshold table
% Print table with the numbers of reliable features according to the
% reliability threshold.

thresholds = [0.5:0.1:0.9, 0.99]';

for tss = 1:nTSS
  % load stats
  if isfile(exp_meta_stats{tss})
    S = load(exp_meta_stats{tss}, 'maxmin_perc', 'mftsLowNaN');
    maxmin_perc{tss} = S.maxmin_perc(S.mftsLowNaN);
  else
    error('File %s is missing!', exp_meta_stats{tss})
  end
end

% calculate numbers of features reliable at least as a defined threshold
tableBase = cellfun(@(perc) ...
                      arrayfun(@(x) sum(~(perc < x)), thresholds), ...
                    maxmin_perc, 'Uni', false);
relTable = table();
for tss = 1:nTSS
  if tss == 1
    relTable.(tssList{tss}) = arrayfun(@(x) ...
                                sprintf('\\textbf{%3d}/%3d', ...
                                  x, numel(maxmin_perc{tss})), ...
                                tableBase{tss}, 'Uni', false);
  else
    relTable.(tssList{tss}) = arrayfun(@(x, y) ...
                                sprintf('\\textbf{%3d}/%3d (%3d/%3d)', ...
                                  y, numel(maxmin_perc{tss}) + numel(maxmin_perc{fullId}), ...
                                  x, numel(maxmin_perc{tss})), ...
                                tableBase{tss}, ...
                                tableBase{tss} + tableBase{fullId}, ...
                                'Uni', false);
  end
end
relTable.Properties.RowNames = arrayfun(@num2str, thresholds, 'Uni', false);

% print table to tex
lt = LatexTable(relTable);
lt.opts.tableColumnAlignment = num2cell('lrrr');
lt.setHeaderRow({'threshold', 'TSS full', 'TSS nearest', 'TSS knn'});
% lt.setColumnFormat({'%s', '%2.2f', []});
lt.opts.tableCaption = [...
  'Proportion of features for individual TSS with robustness greater or equal to the threshold in the first ', ...
  'column. ', ...
  'The proportions in brackets represent $\trainset$-based features ', ...
  'for given TSS (TSS full has only $\archive$-based). ', ...
  'The numbers in bold in the grey row are utilized for the following process.'];
lt.opts.tableLabel = 'featProp';
lt.opts.booktabs = 1;
% add gray background to medoid features
lt.colorizeRowsInGray(5);
lt.toFile(fullfile(tableFolder, 'relTable.tex'));

clear('lt', 'maxmin_perc', 'mftsLowNaN', 'relTable', 'S', 'tableBase', ...
      'tss')

%% Analyse feature reliability
% If the features should be reliable on 90% (99%?) for further
% calculations, each feature should be reliable on
% (0,9)^(1/number of features), i.e., in our case (0,9)^(1/42) - number of
% clusters. Which number of points is the distance between them lower or
% equal to 0.05, i.e. 5% error.

reliability_threshold = 0.05;

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing feature reliability analysis on TSS %s\n'], ...
          fix(clock), tssList{tss})

 if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'observations', ...
             'quantiles_99', ...
             'mfts_names_red');
  else
    error('File %s is missing!', exp_meta_output{tss})
 end

  % set keyword according to TSS method
  if tss == fullId
    keyword = 'archive';
  else
    keyword = 'train';
  end

  % find unique observations values
  obsUni{1} = unique(S.observations(1, :));
  obsUni{2} = unique(S.observations(2, :));

  % use only features with low number of NaNs
  mfts_names_red = S.mfts_names_red;
  
  % init
  nMfts = numel(mfts_names_red);
  % quan_90_diff_break = zeros(nMfts, 1);
  % quan_90_diff_meds = cell(nMfts, 1);
  quan_99_diff_break = zeros(nMfts, 1);
  quan_99_diff_meds = cell(nMfts, 1);
  maxmin_diff_99 = cell(nMfts, 1);
  lowErr_max_int = NaN(nMfts, 2);

  % feature loop
  for ft = 1:nMfts
    close
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Processing feature reliability analysis on TSS %s feature %3d/%d\n'], ...
          fix(clock), tssList{tss}, ft, nMfts)
    % find observation row
    if contains(mfts_names_red{ft}, [keyword, 'test'])
      obsId = 2;
    elseif contains(mfts_names_red{ft}, keyword)
      obsId = 1;
    else
      continue
    end

    % gather quantile differences medians for unique observation values
    for i = 1:numel(obsUni{obsId})
      uniObsId = obsUni{obsId}(i) == S.observations(obsId, :);
      % quan_90_diff_meds{ft}(i) = median(...
      %   S.quantiles_90(ft, uniObsId, 2) - S.quantiles_90(ft, uniObsId, 1));
      quan_99_diff_meds{ft}(i) = median(...
        S.quantiles_99(ft, uniObsId, 2) - S.quantiles_99(ft, uniObsId, 1));
      maxmin_diff_99{ft}(i) = quantile(...
        S.quantiles_99(ft, uniObsId, 2) - S.quantiles_99(ft, uniObsId, 1), 0.99);
    end

    % find the largest interval where the feature error is below 0.05 in
    % 99% cases (quan_99_diff should be identical to max-min)
    lowErr = maxmin_diff_99{ft} < reliability_threshold;
%     if all(lowErr)
%       lowErr_max_int(ft, :) = [obsUni{obsId}(1), obsUni{obsId}(end)];
%     elseif any(lowErr)
%       lowErr_i = find(diff(lowErr));
%       lowErr_n = [lowErr_i, numel(lowErr)] - [0, lowErr_i];
%       [~, lowErr_n_id] = max(lowErr_n(((lowErr(1) == 0)+1) : 2 : end));
%       lowErr_2n_id = lowErr_n_id*2 - (lowErr(1) == 1);
%       if lowErr_2n_id == 1
%         lowErr_max_int(ft, :) = [1, obsUni{obsId}(lowErr_i(lowErr_2n_id))];
%       else
%         lowErr_max_int(ft, :) = [obsUni{obsId}(lowErr_i(lowErr_2n_id - 1) + 1), ...
%                                  obsUni{obsId}(lowErr_i(lowErr_2n_id))];
%       end
%     end
    
%     rel_breakId = 1;
%     rel_break(ft) = obsUni{obsId}(rel_breakId);
%     relObsId = rel_break(ft) <= S.observations(obsId, :);
%     while (sum(S.quantiles_99(ft, relObsId, 2) - S.quantiles_99(ft, relObsId, 1) < reliability_threshold)...
%           / size(S.quantiles_99, 2)) < 0.99 && rel_breakId < numel(obsUni{obsId})
% %       fprintf('Feature %d num of points: %4d\n', ft, rel_break(ft))
%       rel_breakId = rel_breakId + 1;
%       rel_break(ft) = obsUni{obsId}(rel_breakId);
%       relObsId = rel_break(ft) <= S.observations(obsId, :);
%     end

    % plot resulting data
    han = figure();
%     plot(obsUni{obsId}, quan_99_diff_meds{ft})
    plot(obsUni{obsId}, maxmin_diff_99{ft})
    hold on
    % plot(obsUni{obsId}, quan_90_diff_meds{ft})
    % plot reliability threshold
    line([0, obsUni{obsId}(end)], [reliability_threshold, reliability_threshold], 'Color', 'r')
    % plot reliable points
    scatter(obsUni{obsId}(lowErr), reliability_threshold*ones(1, sum(lowErr)), 1, 'g', 'filled')
%     line([lowErr_max_int(1), lowErr_max_int(2)], [reliability_threshold, reliability_threshold], 'Color', [0.4940 0.1840 0.5560])
    xlabel('Number of observations')
    ylabel('max - min')
    title(['TSS ', num2str(tss), ': ', strrep(mfts_names_red{ft}, '_', '\_')])
    legend('max-min Q(0.99)', '0.05 reliability threshold', ...
           'reliable num of points')
%            ['reliable interval [', num2str(lowErr_max_int(1)), ', ', num2str(lowErr_max_int(2)), ']']...
    hold off
    
    % print plot to file
    plotName = fullfile(plotResultsFolder, ['maxmin_tss', num2str(tss), '_ft', num2str(ft)]);
    print2pdf(han, plotName, 1)
%     close
  end
  
  % save results
  save(exp_smsp_quant_anal{tss}, 'reliability_threshold', ...
                                 'lowErr_max_int', ...
                                 'mfts_names_red', 'obsUni')
end

clear('ft', 'i', 'S', 'rel_breakId', ...
      'rel_break', 'relObsId', 'reliability_threshold', ...
      'mfts_names_red', 'obsUni')

close

%% Feature correlation analysis
% Schweizer-Wolff correlation between each pair of available features used
% to calculation of distances between features

nPerms = 5;
threshold = 0.1;

for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing feature correlation analysis on TSS %s\n'], ...
          fix(clock), tssList{tss})

  % load stats
  if isfile(exp_meta_stats{fullId})
    S_stats_full = load(exp_meta_stats{fullId}, 'mftsLowNaN', 'mftsLowErr');
    useFullId = S_stats_full.mftsLowNaN & S_stats_full.mftsLowErr;
  else
    error('File %s is missing!', exp_meta_stats{fullId})
  end
 
  % for each TSS load data for TSS full
  if isfile(exp_meta_output{fullId})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS full output\n'], ...
          fix(clock), tssList{tss})
    S_full = load(exp_meta_output{fullId}, ...
                  'medians', 'inf_mins', 'inf_plus', 'mfts_names_red');
    % use only features with low number of NaNs
    S_full.mfts_names_red = S_full.mfts_names_red(useFullId);
    S_full.medians = S_full.medians(useFullId, :);
    S_full.inf_mins = S_full.inf_mins(useFullId, :);
    S_full.inf_plus = S_full.inf_plus(useFullId, :);
  else
    error('File %s is missing!', exp_meta_output{fullId})
  end
  % load minmax TSS full
  if isfile(exp_meta_minmax{fullId})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS full minmax\n'], ...
          fix(clock), tssList{tss})
    S_mm_full = load(exp_meta_minmax{fullId}, 'ft_min_ninf', 'ft_max_ninf');
    % use only features with low number of NaNs
    S_mm_full.ft_min_ninf = S_mm_full.ft_min_ninf(useFullId);
    S_mm_full.ft_max_ninf = S_mm_full.ft_max_ninf(useFullId);
  else
    error('File %s is missing!', exp_meta_minmax{fullId})
  end

  if tss > 1
      % load stats
      if isfile(exp_meta_stats{tss})
        S_stats_tss = load(exp_meta_stats{tss}, 'mftsLowNaN', 'mftsLowErr');
        useTSSId = S_stats_tss.mftsLowNaN & S_stats_tss.mftsLowErr;
      else
        error('File %s is missing!', exp_meta_stats{tss})
      end

    % load output
    if isfile(exp_meta_output{tss})
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS %s output\n'], ...
          fix(clock), tssList{tss}, tssList{tss})
      S_tss = load(exp_meta_output{tss}, ...
                  'medians', 'inf_mins', 'inf_plus', 'mfts_names_red');
      % use only features with low number of NaNs
      S_tss.mfts_names_red = S_tss.mfts_names_red(useTSSId);
      S_tss.medians = S_tss.medians(useTSSId, :);
      S_tss.inf_mins = S_tss.inf_mins(useTSSId, :);
      S_tss.inf_plus = S_tss.inf_plus(useTSSId, :);
    else
      error('File %s is missing!', exp_meta_output{tss})
    end
    % load minmax
    if isfile(exp_meta_minmax{tss})
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading minmax\n'], ...
          fix(clock), tssList{tss})
      S_mm_tss = load(exp_meta_minmax{tss}, 'ft_min_ninf', 'ft_max_ninf');
      % use only features with low number of NaNs
      S_mm_tss.ft_min_ninf = S_mm_tss.ft_min_ninf(useTSSId);
      S_mm_tss.ft_max_ninf = S_mm_tss.ft_max_ninf(useTSSId);
    else
      error('File %s is missing!', exp_meta_minmax{tss})
    end
  else
    S_tss.medians = [];
    S_tss.inf_mins = [];
    S_tss.inf_plus = [];
    S_tss.mfts_names_red = {};
    S_tss.generations = [];
    S_stats_tss.mftsLowNaN = [];
    S_stats_tss.mftsLowErr = [];
    S_mm_tss.ft_min_ninf = [];
    S_mm_tss.ft_max_ninf = [];
  end

  % concatenate results from TSS full and actual TSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
       'Correlation analysis on TSS %s - cat results with TSS full\n'], ...
      fix(clock), tssList{tss})
  % TODO: precise solution for different numbers of results
  sameId = 1:(min(size(S_full.medians, 2), size(S_tss.medians, 2)));
  medians = [S_full.medians; S_tss.medians(:, sameId)];
  inf_mins = [S_full.inf_mins; S_tss.inf_mins(:, sameId)];
  inf_plus = [S_full.inf_plus; S_tss.inf_plus(:, sameId)];
  ft_min_ninf = [S_mm_full.ft_min_ninf; S_mm_tss.ft_min_ninf];
  ft_max_ninf = [S_mm_full.ft_max_ninf; S_mm_tss.ft_max_ninf];
  mfts_names_red = [S_full.mfts_names_red; S_tss.mfts_names_red];
  mftsLowNaN = [S_stats_full.mftsLowNaN; S_stats_tss.mftsLowNaN];
  mftsLowErr = [S_stats_full.mftsLowErr; S_stats_tss.mftsLowErr];
  nMfts = numel(mfts_names_red);

  medians_minmax = medians;
  % replace NaN's with max or min values, where NaN was placed due to Inf
  % value (should not be necessary when using sigmoid normalization)
  for m = 1:nMfts
    medians_minmax(m, inf_mins(m, :) > 49) = ft_min_ninf(m);
    medians_minmax(m, inf_plus(m, :) > 49) = ft_max_ninf(m);
  end

  % clear all unnecessary variables to free memory
  clear('S_full',  'S_tss', 'S_mm_full', 'S_mm_tss', ...
        'medians', 'inf_mins', 'inf_plus', 'ft_min_ninf', 'ft_max_ninf')

  % calculate Schweizer-Wolff correlations including NaN values pairwise
  % med_corr = corrSchweizer(medians_minmax', 'rows', 'pairwise');
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
       'Correlation analysis on TSS %s - calculating SW correlation\n'], ...
      fix(clock), tssList{tss})
  [med_sim, med_corr, med_2reals, med_2nans] = ...
    nancorr(medians_minmax', 'rows', 'pairwise', 'type', 'Schweizer');

  % create dendrogram from similarity distance (1 - sim) => 0
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
       'Correlation analysis on TSS %s - dendrogram\n'], ...
      fix(clock), tssList{tss})
  % permutation cycle
  for p = 1:nPerms
    med_perm{p} = randperm(nMfts);
    med_link{p} = linkage(squareform(1 - med_sim(med_perm{p}, med_perm{p})));
    med_k_threshold(p) = sum(med_link{p}(:, 3) > threshold);
  end
  % plot dendrogram with red line marking 0.25 distance threshold only for
  % the first permutation
  h_dendr = figure();
  med_dendr = dendrogram(med_link{1}, 0);
  hold on
  line([0, nMfts+1], [threshold, threshold], 'Color', 'r')
  hold off

  % save testing results
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
       'Correlation analysis on TSS %s - saving results\n'], ...
      fix(clock), tssList{tss})
  save(exp_smsp_corr_test{tss}, 'medians_minmax', 'mfts_names_red', ...
                           'mftsLowNaN', 'mftsLowErr', ...
                           'med_corr', 'med_2reals', 'med_2nans', 'med_sim', ...
                           'med_link', 'med_perm', 'med_k_threshold', 'h_dendr')
end

% find number of clusters for next clustering
k = [];
for tss = 1:nTSS
  S = load(exp_smsp_corr_test{tss}, 'med_k_threshold');
  k = [k, S.med_k_threshold];
end
k_res = floor(mean(k));
fprintf('Resulting k = %f', k_res)

% save resulting k to test files
for tss = 1:nTSS
  save(exp_smsp_corr_test{tss}, 'k_res', '-append')
end

clear('m', 'medians', 'medians_minmax', 'nMfts', ...
      'med_sim', 'med_corr', 'med_2reals', 'med_2nans', ...
      'med_link', 'med_dendr', ...
      'S_full',  'S_tss', 'S_mm_full', 'S_mm_tss', 'sameId')

%% Feature correlation clustering
% Cluster features using k-medoids and results from hierarchical clustering

% TSS loop
for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing feature correlation clustering on TSS %s\n'], ...
           fix(clock), tssList{tss})
  if isfile(exp_smsp_corr_test{tss})
    S = load(exp_smsp_corr_test{tss});
  else
    error('File %s is missing!', exp_smsp_corr_test{tss})
  end
  % extract variables
  medians_minmax = S.medians_minmax;
  mfts_names_red = S.mfts_names_red;
  nMfts = numel(mfts_names_red);

  % set number of clusters
  k = S.k_res;
  % set Schweizer-Wolf distance
  corrSWdist = @(x,y) 1-corrSchweizer(x', y', 'rows', 'pairwise');
  % set Schweizer-Wolf distance taking NaNs into account
  corrNanSWdist = @(x,y) 1-nancorr(x', y', 'rows', 'pairwise', 'type', 'Schweizer');

  % % clustering to 60 clusters using k-medoids with wasnan (line 217) set to
  % % false, i.e. NaNs will be removed pairwise
  % fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
  %          'Running k-medoids - NaNs removed pairwise\n'], fix(clock))
  % [corrClusterId_pair, ~, ~, ~, corrMedoidId_pair] = kmedoids2(medians_minmax, k, 'Distance', corrSWdist);

  % % remove NaN columns
  nanCols = any(isnan(medians_minmax), 1);
  % if size(medians_minmax(:, ~nanCols), 2) > 0
  %   % clustering to 60 clusters using k-medoids without columns containing
  %   % at least one NaN
  %   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
  %            'Running k-medoids - without columns with NaN\n'], fix(clock))
  %   [corrClusterId_all, ~, ~, ~, corrMedoidId_all] = kmedoids(medians_minmax(:, ~nanCols), k, 'Distance', corrSWdist);
  % else
  %   fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
  %            'Omitting k-medoids - without columns with NaN ', ...
  %            'because each column contains NaN\n'], fix(clock))
  %   % each column contains NaN
  %   corrClusterId_all = [];
  %   corrMedoidId_all = [];
  % end

  % clustering to 60 clusters using k-medoids with wasnan (line 217) set to
  % false, i.e. NaNs will be taken into account pairwise
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Running k-medoids - NaNs taken into account pairwise\n'], ...
          fix(clock))
  [corrClusterId_nanpair, ~, ~, ~, corrMedoidId_nanpair] = kmedoids2(medians_minmax, k, 'Distance', corrNanSWdist);

  % % count loss of data in 'all' k-medoids strategy
  % dataLoss_all = sum(nanCols)/numel(nanCols);
  % count loss of data in 'pairwise' k-medoids strategy
  dataLoss_pair = sum(sort(sum(isnan(medians_minmax), 2))'.*(0:(nMfts-1))) / nchoosek(nMfts,2) / numel(nanCols);

  % identify selected medoids in mfts list
  mftsIds_red = find(mftsIds);
  % use only features with low number of NaNs
  mftsIds_red = mftsIds_red(logical(S.mftsLowNaN & S.mftsLowErr));
  mftsMedoidIds = mftsIds_red(corrMedoidId_nanpair);

  save(exp_smsp_corr_cluster{tss}, 'k', 'corrSWdist', 'corrNanSWdist', ...
                              'corrClusterId_nanpair', 'corrMedoidId_nanpair', ...
                              'dataLoss_pair', ...
                              'mfts_names_red', 'mftsMedoidIds', ...
                              'lastArchFtId')
  %                             'corrClusterId_pair', 'corrMedoidId_pair', ...
  %                             'corrClusterId_all', 'corrMedoidId_all', ...
  %                             'dataLoss_all',
end

clear('corrSWdist', 'corrNaNSWdist', ...
      'corrClusterId_pair', 'corrMedoidId_pair', ...
      'corrClusterId_all', 'corrMedoidId_all', ...
      'corrClusterId_nanpair', 'corrMedoidId_nanpair', ...
      'dataLoss_all', 'dataLoss_pair', 'k', 'nanCols', 'S', 'tss', ...
      'mfts_names_red')

%% Feature correlation clustering analysis
% Analyse results of feature clustering

% maximal number of features in table per page
maxFtsPerPage = 40;

% TSS loop
for tss = 1:nTSS
  if isfile(exp_smsp_corr_cluster{tss})
    S_corr = load(exp_smsp_corr_cluster{tss});
  else
    error('File %s is missing!', exp_smsp_corr_cluster{tss})
  end

  % load TSS full simple stats
  if isfile(exp_meta_stats{fullId})
    S_stats_full = load(exp_meta_stats{fullId}, ...
                        'mftsLowNaN', 'mftsLowErr', 'maxmin_perc_rel', ...
                        'maxmin_perc', 'mfts_names_red');
    useFullId = S_stats_full.mftsLowNaN & S_stats_full.mftsLowErr;
  else
    error('File %s is missing!', exp_meta_stats{fullId})
  end
  % add actual TSS results
  if tss ~= fullId
    if isfile(exp_meta_stats{tss})
      S_stats_tss = load(exp_meta_stats{tss}, ...
                        'mftsLowNaN', 'mftsLowErr', 'maxmin_perc_rel', ...
                        'maxmin_perc', 'mfts_names_red');
      useTSSId = S_stats_tss.mftsLowNaN & S_stats_tss.mftsLowErr;
      % cat TSS full and actual TSS results
      maxmin_perc = [S_stats_full.maxmin_perc(:); S_stats_tss.maxmin_perc(:)];
      maxmin_perc_rel = [S_stats_full.maxmin_perc_rel(:); S_stats_tss.maxmin_perc_rel(:)];
      mfts_names_red = [S_stats_full.mfts_names_red; S_stats_tss.mfts_names_red];
    else
      error('File %s is missing!', exp_meta_stats{tss})
    end
  else
    maxmin_perc = S_stats_full.maxmin_perc(:);
    maxmin_perc_rel = S_stats_full.maxmin_perc_rel(:);
    mfts_names_red = S_stats_full.mfts_names_red;
    useTSSId = logical([]);
  end

  % load TSS full NaN analysis
  if isfile(exp_smsp_nan_anal{fullId})
    S_nan_full = load(exp_smsp_nan_anal{fullId});
  else
    error('File %s is missing!', exp_smsp_nan_anal{fullId})
  end
  % add actual TSS results
  if tss ~= fullId
    if isfile(exp_smsp_nan_anal{tss})
      S_nan_tss = load(exp_smsp_nan_anal{tss});
      % cat TSS full and actual TSS results
      nan_break = [S_nan_full.nan_break; ...
                   S_nan_tss.nan_break];
      nan_break_rel = [S_nan_full.nan_break(useFullId); ...
                       S_nan_tss.nan_break(useTSSId)];
    else
      error('File %s is missing!', exp_smsp_nan_anal{tss})
    end
  else
    nan_break = S_nan_full.nan_break;
    nan_break_rel = S_nan_full.nan_break(useFullId);
  end
  
  % TODO: dimension test was probably performed without TSS full and actual
  % TSS unification => check and perhaps unify here
  % TSS full dimension test result
  if isfile(exp_smsp_dimension_test{fullId})
    S_dim_full = load(exp_smsp_dimension_test{fullId}, ...
                      'dims', ...
                      'wilcoxon_uniq_combs', 'wilcoxon_uniq_combs_id');
  else
    error('File %s is missing!', exp_smsp_dimension_test{fullId})
  end
  % add actual TSS results
  if tss ~= fullId
    if isfile(exp_smsp_dimension_test{tss})
      S_dim_tss = load(exp_smsp_dimension_test{tss}, ...
                      'wilcoxon_uniq_combs', 'wilcoxon_uniq_combs_id');
    else
      error('File %s is missing!', exp_smsp_dimension_test{tss})
    end
  else
    S_dim_tss.wilcoxon_uniq_combs = [];
    S_dim_tss.wilcoxon_uniq_combs_id = [];
  end
  % prepare necessary variables and concatenate actual TSS and TSS full
  % dimensional testing results
  dims = S_dim_full.dims;
  nMfts = numel(mfts_names_red);
  nMfts_rel = numel(S_corr.mfts_names_red);
  
  % dimension test results
  dimCombs = nchoosek(dims, 2);
  dimCombsStr = arrayfun(@(x) ...
                  sprintf('(%d,\\,%d)', dimCombs(x, 1), dimCombs(x, 2)), ...
                  1:size(dimCombs, 1), 'uni', false);
  % TSS full combinations
  dimCombsAccepted = arrayfun(@(x) ...
    strjoin( ...
      dimCombsStr( ...
        S_dim_full.wilcoxon_uniq_combs( ...
          S_dim_full.wilcoxon_uniq_combs_id(x), ...
          :) ...
        ), ...
      ', ' ...
    ), ...
    1:numel(S_dim_full.wilcoxon_uniq_combs_id), ...
    'UniformOutput', false)';
  % add actual TSS combinations
  dimCombsAccepted = [dimCombsAccepted; ...
    arrayfun(@(x) ...
      strjoin( ...
        dimCombsStr( ...
          S_dim_tss.wilcoxon_uniq_combs( ...
            S_dim_tss.wilcoxon_uniq_combs_id(x), ...
            :) ...
          ), ...
        ', ' ...
      ), ...
      1:numel(S_dim_tss.wilcoxon_uniq_combs_id), ...
      'UniformOutput', false)' ...
    ];

  % table mfts notation
  mftsSplit = cellfun(@(x) strsplit(x, '_'), mfts_names_red, 'Uni', false);
  setColumn   = cell(nMfts, 1);
  transColumn = cell(nMfts, 1);
  classColumn = cell(nMfts, 1);
  classStarts = [];
  for m = 1:nMfts
    % set notation
    if any(strcmp(mftsSplit{m}{1}, ...
        {'archive', 'archivetest', 'traintest', 'train'}))
      switch mftsSplit{m}{1}
        case 'archive'
          setColumn{m} = '\archive';
        case 'archivetest'
          setColumn{m} = '\archivepred';
        case 'train'
          setColumn{m} = '\trainset';
        case 'traintest'
          setColumn{m} = '\trainpredset';
      end
      % remove set notation
      mftsSplit{m}(1) = [];
    end
    % transformation notation
    if any(strcmp(mftsSplit{m}{1}, ...
        {'none', 'cma'}))
      switch mftsSplit{m}{1}
        case 'none'
          transColumn{m} = '\transnone';
        case 'cma'
          transColumn{m} = '\transcma';
      end
      % remove transformation notation
      mftsSplit{m}(1) = [];
    end
    % feature class notation
    if any(strcmp(mftsSplit{m}{1}, ...
        {'cmaes', 'dispersion', 'ela', 'infocontent', 'nearest'}))
       switch mftsSplit{m}{1}
         case 'cmaes'
           classColumn{m} = '\text{\ftCMA}';
         case 'dispersion'
           classColumn{m} = '\text{\ftDisp}';
         case 'ela'
           switch mftsSplit{m}{2}
             case 'distribution'
               classColumn{m} = '\text{\ftyDis}';
             case 'levelset'
               classColumn{m} = '\text{\ftLevel}';
             case 'metamodel'
               classColumn{m} = '\text{\ftMM}';
           end
           % remove ela set notation
           mftsSplit{m}(1) = [];
         case 'infocontent'
           classColumn{m} = '\text{\ftInfo}';
         case 'nearest'
           classColumn{m} = '\text{\ftNBC}';
           % remove 'nearest' part of set notation
           mftsSplit{m}(1) = [];
       end
      % remove set notation
      mftsSplit{m}(1) = [];
    end
    % unite the rest
    mftsSplit{m} = strjoin(mftsSplit{m}, '\\_');
    % mark different classes
    if m > 1 && ~isequal(classColumn{m}, classColumn{m-1})
      classStarts(end+1) = m;
    end
  end
  tableMftsNotation = cellfun(...
    @(x, y, z, w) sprintf('$\\tableFeat{%s}{%s}{%s}{%s}$', x, y, z, w), ...
    mftsSplit, setColumn, transColumn, classColumn, 'Uni', false);

%   % all feature properties page cycle
%   pages = ceil(nMfts/maxFtsPerPage);
%   useInTSS = [useFullId; useTSSId];
%   for p = 1:pages
%     rowIds = (p-1)*maxFtsPerPage+1 : min(p*maxFtsPerPage, nMfts);
%     % table with correlation clustering ids
%     featurePropTable = table(...
%                              nan_break(rowIds), ...
%                              100*maxmin_perc(rowIds), ...
%                              dimCombsAccepted(rowIds), ...
%                              'RowNames', tableMftsNotation(rowIds));
% 
%     % print table to tex
%     lt = LatexTable(featurePropTable);
%     lt.opts.tableColumnAlignment = num2cell('lrrl');
%     lt.setHeaderRow({'', '$N_\nanout$', 'rob.(\%)', '$(\dm_i, \dm_j)$'});
%     lt.setColumnFormat({'%d', '%2.2f', []});
%     lt.opts.midLines = [1, ... % first line
%                         classStarts(classStarts < max(rowIds) & ...
%                                     classStarts > min(rowIds)) ...
%                           - min(rowIds) + 1];
%     lt.opts.tableCaption = sprintf([...
%       'TSS %s features (%d/%d). ', ...
%       'Features are grouped according to their feature sets (separated by horizontal lines). ', ...
%       'Features with less than 25\\%% of values equal to $\\nanout$ and robustness greater than 0.9, are marked as gray lines. ', ...
%       '$N_\\nanout$ denotes the lowest measured number of points from which at most 1\\%% of feature calculations resulted in $\\nanout$. ', ...
%       'The $(\\dm_i, \\dm_j)$ column shows the pairs of feature dimensions ', ...
%       'for which the two-sided Wilcoxon signed rank test with the Bonferroni-Holm correction ', ...
%       'does not reject the hypothesis of equality of median feature values, ', ...
%       'at the family-wise level 0.05 for each individual feature.'...
%       ], tssList{tss}, p, pages);
%     lt.opts.tableLabel = sprintf('featProp_%s_%d', tssList{tss}, p);
%     lt.opts.booktabs = 1;
%     % add gray background to medoid features
%     lt.colorizeRowsInGray(useInTSS(rowIds));
%     lt.toFile(sprintf('%s_%d.tex', exp_smsp_feat_table{tss}, p));
%   end
  
  % table to appendix - contains multicolumns with individual sample sets
  decomposedNotation = extractBetween(tableMftsNotation, '{', '}');
  appNotId = strcmp(decomposedNotation(:, 2), '\archive') & ...
             (strcmp(decomposedNotation(:, 3), '\transnone') | ...
              strcmp(decomposedNotation(:, 3), '')) | ...
             cellfun(@isempty, decomposedNotation(:, 2));
  % rename dimension and observations
  tableMftsNotation(strcmp(decomposedNotation(:, 1), 'dimension'), 1) = ...
    strrep(tableMftsNotation(strcmp(decomposedNotation(:, 1), 'dimension'), 1), ...
           'dimension', 'dim');
  tableMftsNotation(strcmp(decomposedNotation(:, 1), 'observations'), 1) = ...
    strrep(tableMftsNotation(strcmp(decomposedNotation(:, 1), 'observations'), 1), ...
           'observations', 'obs');
%   tableFeat(strcmp(decomposedNotation, 'observations'), 1) = {'obs'};
  % remove sample set (archive) notation
  appendixNotation = strrep(tableMftsNotation(appNotId), '\archive', '');
  nAppNotation = numel(appendixNotation);
  % gain sample sets without transformation (no sample set excluded)
%   sampleSet = unique(decomposedNotation(:, 2));
%   sampleSet = sampleSet(2:end);
  if tss == fullId
    sampleSet = {'\archive', '\archivepred'};
    tableSampleSet = {'\archive', '\trans{\archive}', ...
                      '\archivepred', '\trans{\archivepred}'};
  else
    sampleSet = {...'\archive', '\archivepred', ...
                 '\trainset', '\trainpredset'};
    tableSampleSet = {...'\archive', '\trans{\archive}', ...
                      ...'\archivepred', '\trans{\archivepred}', ...
                      '\trainset', '\trans{\trainset}', ...
                      '\trainpredset', '\trans{\trainpredset}'};
  end
  nSampleSet = numel(sampleSet);
  % gain transformation settings (empty transformation excluded)
%   transSet = unique(decomposedNotation(:, 3));
%   transSet = transSet(2:end);
  transSet = {'\transnone', '\transcma'};
  nTransSet = numel(transSet);
  % usage of individual features
  useInTSS = [useFullId; useTSSId];
  % init
  nan_break_app = NaN(nAppNotation, nSampleSet*nTransSet);
  maxmin_perc_app = NaN(nAppNotation, nSampleSet*nTransSet);
  dimCombsAccepted_app = {};
  dimCombsAccepted_app(1:nAppNotation, 1:nSampleSet*nTransSet) = {NaN};
  useInTSS_app = false(nAppNotation, nSampleSet*nTransSet);
  % column cycle
  for s = 1:nSampleSet
    for t = 1:nTransSet
      sampleSetId = strcmp(sampleSet{s}, decomposedNotation(:, 2)) & ...
                    (strcmp(transSet{t}, decomposedNotation(:, 3)) | ...
                     strcmp('', decomposedNotation(:, 3)));
      rowPositionId = ismember(decomposedNotation(appNotId, 1), ...
                               decomposedNotation(sampleSetId, 1));
      nan_break_app(rowPositionId, (s-1)*nTransSet + t) = nan_break(sampleSetId);
      maxmin_perc_app(rowPositionId, (s-1)*nTransSet + t) = maxmin_perc(sampleSetId);
      dimCombsAccepted_app(rowPositionId, (s-1)*nTransSet + t) = dimCombsAccepted(sampleSetId);
      useInTSS_app(rowPositionId, (s-1)*nTransSet + t) = useInTSS(sampleSetId);
    end
  end
  % fill in independent features
  indepId = cellfun(@isempty, decomposedNotation(:, 2));
  rowIndepId = ismember(decomposedNotation(appNotId, 1), ...
                        decomposedNotation(indepId, 1));
  nan_break_app(rowIndepId, : ) = repmat(nan_break(indepId), 1, nSampleSet*nTransSet);
  maxmin_perc_app(rowIndepId, : ) = repmat(maxmin_perc(indepId), 1, nSampleSet*nTransSet);
  dimCombsAccepted_app(rowIndepId, : ) = repmat(dimCombsAccepted(indepId), 1, nSampleSet*nTransSet);
  useInTSS_app(rowIndepId, :) = repmat(useInTSS(indepId), 1, nSampleSet*nTransSet);
  
  % all feature properties page cycle
  maxFtsPerAppTable = 64;
  pages = ceil(nAppNotation/maxFtsPerAppTable);
  maxSetsPerTable = 4;
  columnTables = ceil(nSampleSet*nTransSet/maxSetsPerTable);
  for c = 1:columnTables
    for p = 1:pages
      rowIds = (p-1)*maxFtsPerAppTable+1 : min(p*maxFtsPerAppTable, nAppNotation);
      % table content
      tabContent = {};
      grayContent = useInTSS_app(rowIds, (c-1)*maxSetsPerTable + 1 : c*maxSetsPerTable);
      for s = (c-1)*maxSetsPerTable + 1 : c*maxSetsPerTable
        tabContent = [tabContent, {...
                       nan_break_app(rowIds, s), ...
                       100*maxmin_perc_app(rowIds, s), ...
                       dimCombsAccepted_app(rowIds, s) ...
                     }];
      end
      % table with correlation clustering ids
      featurePropTable = table(...
                               tabContent{:}, ...
                               'RowNames', appendixNotation(rowIds));

      % print table to tex
      lt = LatexTable(featurePropTable);
      lt.opts.tableColumnAlignment = [{'L{25mm}'}, repmat(num2cell('R{5.5mm}R{6.5mm}L{7mm}'), 1, maxSetsPerTable)];
      lt.setHeaderRow([{''}, repmat({'$N_\nanout$', 'rob.(\%)', '$(\dm_i, \dm_j)$'}, 1, maxSetsPerTable)]);
      lt.setColumnFormat(repmat({'%d', '%2.2f', []}, 1, maxSetsPerTable));
      lt.opts.midLines = [1, ... % first line
                          classStarts(classStarts < max(rowIds) & ...
                                      classStarts > min(rowIds)) ...
                            - min(rowIds) + 1];
      % if pages*columnTables > 1
      %   tableCounter = sprintf(' (%d/%d)', (c-1)*pages + p, pages*columnTables);
      % else
      %   tableCounter = '';
      % end
      if tss == fullId
        tableCounter = ' ($\archive=\trainset$ for TSS full)';
      else
        tableCounter = [' (only $\trainset$-based and sample set indepent,', ...
                        ' $\archive$-based are identical to TSS full', ...
                        ' in Table~\ref{table:featProp_full_1})'];
      end
      lt.opts.tableCaption = sprintf([...
        'TSS %s features%s. ', ...
        'Features are grouped according to their feature sets (separated by horizontal lines). ', ...
        'Features with less than 25\\%% of values equal to $\\nanout$ and robustness greater than 0.9, are in gray. ', ...
        '$N_\\nanout$ denotes the lowest measured number of points from which at most 1\\%% of feature calculations resulted in $\\nanout$ ', ...
        '($N_\\nanout = 0$ for sample set independent $\\feat{}$). ', ...
        'The $(\\dm_i, \\dm_j)$ column shows the pairs of feature dimensions ', ...
        'for which the two-sided Wilcoxon signed rank test with the Bonferroni-Holm correction ', ...
        'does not reject the hypothesis of equality of median feature values, ', ...
        'at the family-wise level 0.05 for each individual feature.'...
        ], tssList{tss}, tableCounter);
      lt.opts.tableLabel = sprintf('featProp_%s_%d', tssList{tss}, (c-1)*pages + p);
      lt.opts.booktabs = 1;
      % add gray background to medoid features
%       lt.colorizeRowsInGray(useInTSS(rowIds));
%       lt.colorizeSubMatrixInGray(values, row, col, minGray, maxGray);
      for rc = 1:size(grayContent, 1)
        for cc = 1:size(grayContent, 2)
          if grayContent(rc, cc)
            lt.colorizeSubMatrixInGray([1, 1, 1], rc, 3*(cc-1)+1, 0.85, 0.85);
          end
        end
      end
      % add sample set header row
      ssHeaderRow = strjoin(cellfun(@(x) ['\multicolumn{3}{c}{$', x, '$}'], ...
        tableSampleSet((c-1)*maxSetsPerTable + 1 : c*maxSetsPerTable), ...
        'Uni', false), ' & ');
      latexRows = lt.toStringRows;
      latexRows = [{'\newcommand{\PreserveBackslash}[1]{\let\temp=\\#1\let\\=\temp}', ...
                    '\newcolumntype{C}[1]{>{\PreserveBackslash\centering}p{#1}}', ...
                    '\newcolumntype{R}[1]{>{\PreserveBackslash\raggedleft}p{#1}}', ...
                    '\newcolumntype{L}[1]{>{\PreserveBackslash\raggedright}p{#1}}'}';
                   latexRows(1:4); ...
                   {'\resizebox{\textwidth}{!}{%'}; ...
                   latexRows(5:6); ...
                   {['{} & ', ssHeaderRow,' \\']}; ...
                   arrayfun(@(x, y) sprintf('\\cmidrule(lr){%d-%d}', x, y), ...
                   2:3:3*maxSetsPerTable, (2:3:3*maxSetsPerTable)+2, 'Uni', false)'; ...
                   latexRows(7:end-1); ...
                   {'}'}; ...
                   latexRows(end)];
      % save the result in the file
      fid = fopen(sprintf('%s_%d.tex', exp_smsp_feat_table{tss}, (c-1)*pages + p), 'w');
      for i = 1:length(latexRows)
        fprintf(fid, '%s\n', latexRows{i});
      end
      fclose(fid);
%       lt.toFile(sprintf('%s_%d.tex', exp_smsp_feat_table{tss}, (c-1)*pages + p));
    end
  end
  
  % reduce features to clustered features
  tableMftsNotation_red = tableMftsNotation([useFullId; useTSSId]);
  
  % sort metafeatures to clusters
%   [~, ia, ic] = unique(S_corr.corrClusterId_nanpair);
  [~, sortClId] = sort(S_corr.corrClusterId_nanpair);

  % is medoid id
  isMed = ismember((1:nMfts_rel)', S_corr.corrMedoidId_nanpair);

  % sort table variables
  nan_break_rel = nan_break_rel(sortClId);
  maxmin_perc_rel = 100*maxmin_perc_rel(sortClId);
  tableMftsNotationSorted{tss} = tableMftsNotation_red(sortClId);
  corrClusterId_nanpairSorted{tss} = S_corr.corrClusterId_nanpair(sortClId);
  isMedSorted{tss} = isMed(sortClId);

  % page cycle
  switch tss
    case fullId
      rowIdsTss = {[30:33, 4:17], [18:29, 1:3]};
    case nearId
      rowIdsTss = {2:25, [26:45, 1], [59:60, 46:50], 51:58};
    case knnId
      rowIdsTss = {[1:24, 53:55, 27:30], [31:52, 25:26, 56:59]};
  end
%   if tss == nearId
%     ftsPerColumn = [25, 20, 7, 8];
%   else
%     ftsPerColumn = [ceil(nMfts_rel/2), ...
%                    floor(nMfts_rel/2)];
%   end
  columns = numel(rowIdsTss);
  for c = 1:columns
    rowIds = rowIdsTss{c};
%     if c == 1
%       rowIds = 1:ftsPerColumn(1);
%     else
%       rowIds = sum(ftsPerColumn(1:c-1))+1 : sum(ftsPerColumn(1:c));
%     end
    % table with correlation clustering ids
    corrClusterTable = table(...
                             nan_break_rel(rowIds), ...
                             maxmin_perc_rel(rowIds), ...
                             'RowNames', tableMftsNotationSorted{tss}(rowIds));

    % print table to tex
    lt = LatexTable(corrClusterTable);
    lt.opts.tableColumnAlignment = num2cell('lrr');
    lt.setHeaderRow({'', '$N_\nanout$', 'rob.(\%)'});
    lt.setColumnFormat({'%d', '%2.2f'});
    [~, lt.opts.midLines] = unique(corrClusterId_nanpairSorted{tss}(rowIds));
    lt.opts.booktabs = 1;
    lt.opts.latexHeader = false;
    % add gray background to medoid features
    lt.colorizeRowsInGray(isMedSorted{tss}(rowIds));
    lt.toFile(sprintf('%s_%d.tex', exp_smsp_corr_dim_table{tss}, c));
  end

  % save mfts notation
  save(exp_smsp_corr_cluster{tss}, ...
         'tableMftsNotation', 'tableMftsNotation_red', '-append')
end

% cluster table - only clustering
nMaxFeat = max(cellfun(@numel, isMedSorted));
clusterNumberCol = cell(1, 3);
clStarts = cell(1, 3);
% cluster numbering cycle
for tss = 1:nTSS
  clusterNumberCol{tss} = cell(nMaxFeat, 1);
  [clUni, clStarts{tss}] = unique(corrClusterId_nanpairSorted{tss});
  clusterNumberCol{tss}(clStarts{tss}) = ...
    arrayfun(@(x, y) ...
      sprintf( ...
        '\\multirow{%d}{*}{%d}', ...
        x, ...
        y ...
      ), ...
      diff([clStarts{tss}; numel(corrClusterId_nanpairSorted{tss}) + 1]), ...
      clUni, ...
      'Uni', false);
end
clusterTable = table(...
                     clusterNumberCol{fullId}, ...
                     [tableMftsNotationSorted{1}; cell(27, 1)], ...
                     clusterNumberCol{nearId}, ...
                     tableMftsNotationSorted{2}, ...
                     clusterNumberCol{knnId}, ...
                     [tableMftsNotationSorted{3}; cell(1)] ...
                    );
% clusterTable = table(...
%                      [arrayfun(@num2str, corrClusterId_nanpairSorted{fullId}, 'Uni', 0); cell(27, 1)], ...
%                      [tableMftsNotationSorted{1}; cell(27, 1)], ...
%                      arrayfun(@num2str, corrClusterId_nanpairSorted{nearId}, 'Uni', 0), ...
%                      tableMftsNotationSorted{2}, ...
%                      [arrayfun(@num2str, corrClusterId_nanpairSorted{knnId}, 'Uni', 0); cell(1, 1)], ...
%                      [tableMftsNotationSorted{3}; cell(1)] ...
%                     );
% print table to tex
lt = LatexTable(clusterTable);
lt.opts.tableColumnAlignment = num2cell('rlrlrl');
lt.setHeaderRow({'\multicolumn{2}{c}{TSS full}', '', ...
                 '\multicolumn{2}{c}{TSS nearest}', '', ...
                 '\multicolumn{2}{c}{TSS knn}', ''});
% [~, lt.opts.midLines] = unique(corrClusterId_nanpairSorted{tss});
lt.opts.booktabs = 1;
lt.opts.latexHeader = false;
% add gray background to medoid features
for tss = 1:nTSS
  arrayfun(@(x) ...
    lt.prependFormatXY(x, 2*tss, '\\cellcolor[gray]{0.85}'), find(isMedSorted{tss})...
  )
end
latexRows = cell(2*numel(lt.toStringRows), 1);
cmidRows = cell(size(lt.toStringRows));
% create cmidrules
for tss = 1:nTSS
  cmidRows(clStarts{tss}(2:end) + 3) = cellfun(@(x) [x, ...
                                         sprintf('\\cmidrule(lr){%d-%d}', ...
                                                 2*(tss-1) + 1, 2*(tss-1) + 2)], ...
                                       cmidRows(clStarts{tss}(2:end) + 3), 'Uni', false);
end
% concatenate table rows and mid rows
latexRows(1:2:end) = lt.toStringRows;
latexRows(2:2:end) = cmidRows;
lt.toFile(fullfile(tableFolder, 'clusterTable.tex'));
% save the result to the file
fid = fopen(fullfile(tableFolder, 'clusterTable.tex'), 'w');
for i = 1:length(latexRows)
  fprintf(fid, '%s\n', latexRows{i});
end
fclose(fid);

clear('c', 'm', 'nMfts', 'S_dim_full', 'S_dim_tss')

%% Feature means and variances
% The following graphs show dependencies of feature means and variances on
% the number of observations and data densities.
% Means, variances, number of NaNs, and Infs were calculated using
% features
% computed on 100 
% sampled datasets utilizing the same CMA-ES distribution (100 sampled
% _archives_, _training_, and _testing sets_). These sample statistics were
% grouped according to quantiles of number of observations and distribution
% density. Each graph shows dependency of 0.05, 0.5, and 0.95 quantile of 
% qroup feature mean or variance (blue color, right axis) 
% on 1000 quantiles of data density or
% number of observation. In case of non-zero number of NaNs or Infs, the
% group median number of NaNs or Infs (red or green color, left axis) 
% dependencies are also shown.
% To identify possible influence of dimension values in idividual groups,
% the presence of data with particular dimension in quantile group is 
% marked in the graph background (violet color, left axis).
%
% Feature names follow the pattern *{type of data}_{feature
% class}_{feature}*, where *{type of data}* refers to the set of data used
% for feature computation (_archive_, _training_, or _training + testing 
% set_),
% *{feature class}* denotes the class of features (_CMA-ES_ based, 
% _Dispersion_, _ELA y-Distribution_, _ELA Levelset_, _ELA Metamodel_, 
% _Information
% Content_, and _Nearest Better Clustering_), and *{feature}* denotes the
% specific feature.

sampleSetId = zeros(numel(mfts_names), 1);
for o = 1:size(observations, 1)
  sampleSetId(contains(mfts_names, sampleSets{o})) = o;
end

gen_mfts_names_prt = mfts_names_prt(sampleSetId == 0);

% general features
% for ft = 1 :sum(sampleSetId == 0)
%   %%
%   close all
%   
%   % general feature
%   for o = 1:size(observations, 1)
%     figure()
%     % plot means on obs
%     plot(un_observations{o}, general_obs_means(ft + (o-1)*nGenFeat, 1:length(un_observations{o})))
% 
%     hold on
%     title([sampleSets{o}(1:end-1), '\_', gen_mfts_names_prt{ft}])
% %     legend('Values')
%     xlabel('Number of observations')
%     ylabel('Mean')
%     hold off
%   end
%   
% end

for ft = 1:numel(mfts_names)
  %%
  close all
  
  % nan, Inf, -Inf medians
  med_obs_nans = obs_nans(ft, :, 2);
  med_obs_inf_plus = obs_inf_plus(ft, :, 2);
  med_obs_inf_mins = obs_inf_mins(ft, :, 2);
  
  med_dens_nans = dens_nans(ft, :, 2);
  med_dens_inf_plus = dens_inf_plus(ft, :, 2);
  med_dens_inf_mins = dens_inf_mins(ft, :, 2);
  
  if sampleSetId(ft) > 0
    
    fprintf('%s\n', mfts_names{ft})
    
    % observations
    mftPropPlot(shiftdim(obs_means(ft, :, :), 1), med_obs_nans, ...
                med_obs_inf_plus, med_obs_inf_mins, ...
                'Dimensions', shiftdim(obs_dims(sampleSetId(ft), :, :), 1), ...
                'DimTypes', dims, ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 mean', '0.5 mean', '0.95 mean'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', obs_obs(sampleSetId(ft), :), ...
                'XLabelName', 'number of observations');
              
    mftPropPlot(shiftdim(obs_vars(ft, :, :), 1), med_obs_nans, ...
                med_obs_inf_plus, med_obs_inf_mins, ...
                'Dimensions', shiftdim(obs_dims(sampleSetId(ft), :, :), 1), ...
                'DimTypes', dims, ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 var', '0.5 var', '0.95 var'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', obs_obs(sampleSetId(ft), :), ...
                'XLabelName', 'number of observations');
              
    % densities
    mftPropPlot(shiftdim(dens_means(ft, :, :), 1), med_dens_nans, ...
                med_dens_inf_plus, med_dens_inf_mins, ...
                'Dimensions', shiftdim(dens_dims(sampleSetId(ft), :, :), 1), ...
                'DimTypes', dims, ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 mean', '0.5 mean', '0.95 mean'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', dens_dens(sampleSetId(ft), :), ...
                'XLabelName', 'density of observations');
              
    mftPropPlot(shiftdim(dens_vars(ft, :, :), 1), med_dens_nans, ...
                med_dens_inf_plus, med_dens_inf_mins, ...
                'Dimensions', shiftdim(dens_dims(sampleSetId(ft), :, :), 1), ...
                'DimTypes', dims, ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 var', '0.5 var', '0.95 var'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', dens_dens(sampleSetId(ft), :), ...
                'XLabelName', 'density of observations');

  end
end

%% 
 
% % Distribution plots
% 
% if isfile(exp_meta_inf) && false
% 
%   load(exp_meta_inf)
%   
%   for ft = 1:numel(mfts_names)
%     %%
%     close all
%     % plot quantiles
%     figure()
%     plot(q_levels, q_ft(ft, :))
% 
%     hold on
%     % inf values
%     if inf_plus(ft) > 0
%       ylim = get(gca, 'YLim');
%       q_inf = q_ft(ft,:) == Inf;
%       scatter(q_levels(q_inf), ylim(2)*ones(1, sum(q_inf)), 'r', 'filled')
%     elseif inf_min(ft) > 0
%       ylim = get(gca, 'YLim');
%       q_inf = q_ft(ft,:) == -Inf;
%       scatter(q_levels(q_inf), ylim(1)*ones(1, sum(q_inf)), 'r', 'filled')
%     end
%     title(mfts_names_prt{ft})
%     if any(isinf(q_ft(ft, :)))
%       legend('Values', 'Infs')
%     else
%       legend('Values')
%     end
%     xlabel('Quantile')
%     ylabel('Value')
% 
%     hold off
% 
%   end
% else
%   warning('Cannot plot distribution. File %s does not exist', exp_meta_inf)
% end

%% 

% Graphical view

% for ft = 1:numel(mfts_names)
%   
%   %%
%   
%   close all
%   
%   fprintf('\n')
%   
%   nNaNId = ~isnan(obs_vars(ft, :));
%   plot_obs = un_observations(nNaNId);
%  
%   figure()
%   plot(plot_obs, obs_means(ft, nNaNId))
%   title(mfts_names_prt{ft})
%   xlabel('Number of observations')
%   ylabel('Mean')
%   
%   figure()
%   plot(plot_obs, obs_vars(ft, nNaNId))
%   title(mfts_names_prt{ft})
%   xlabel('Number of observations')
%   ylabel('Variance')
%   
% end

%% 

% Visual inspection
% Medians (thick lines) and quartiles (thin and dash-dotted lines) of GP 
% model RDE dependency on individual metafeatures.

% if printScriptMess
%   fprintf('Starting visual inspection\n')
% end
% nPointsToPlot = 200;
% logBound = 5;
% medianLineWidth = 1.8;
% quartileLineWidth = 1;
% medianLineStyle = '-';
% quartileLineStyle = '-.';
% 
% kerColor = getAlgColors([1, 2, 3, 12, 10, 11, 8, 5]) / 255;
% 
% nExtF = 55;
% mfts_order = [14, 16:18, 20:21, ... identical
%               15, 80, 138, ... observations
%               19, 78, 136, ... cma_mean_dist
%               22, 79, 137, ... cma_lik
%               reshape(22 + repmat([0,58,116]', 1, nExtF) ...
%                 + [(1:nExtF); (1:nExtF); (1:nExtF)], 1, 3*nExtF)];
% 
% close all
%  
% % metafeaturePlot(full_mfts_vis(:, mfts_order), full_mfts_vis(:, 6:13), ...
% %                 'DataColor', kerColor, ...
% %                 'DataNames', modelLabels, ...
% %                 'MftsNames', full_mfts_vis.Properties.VariableNames(mfts_order), ...
% %                 'NValues', 200, ...
% %                 'LogBound', 2, ...
% %                 'QuantileRange', [0.05, 0.95], ...
% %                 'MedianLW', 1.8, ...
% %                 'QuartileLW', 1, ...
% %                 'MedianLS', '-', ...
% %                 'QuartileLS', '-.' ...
% %   );
% 
% % print archive_cma_lik
% mftId = 22;
% mftName = ' ';
% pdfNames = {fullfile(plotResultsFolder,     'archive_cma_lik.pdf')};
% 
% han = metafeaturePlot(table2array(full_mfts_vis(:, 22)), full_mfts_vis(:, 6:13), ...
%                 'DataColor', kerColor, ...
%                 'DataNames', modelLabels, ...
%                 'MftsNames', mftName, ...
%                 'NValues', 200, ...
%                 'LogBound', 2, ...
%                 'QuantileRange', [0.05, 0.95], ...
%                 'MedianLW', 1.8, ...
%                 'QuartileLW', 1, ...
%                 'MedianLS', '-', ...
%                 'QuartileLS', '-.' ...
%   );
% 
% print2pdf(han, pdfNames, 1)

%% 

% Correlation analysis
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
% Distribution plots

% close all
% q_bound = [0.05, 0.95];
% 
% for r = 2 % 2:4
%   for c = 20% 1:53 + 5*sign(4-r)
%     distr_all = ks_res.Distributions{r,c};
%     q_d_all = quantile(distr_all, q_bound);
%     
%     d_all_show = distr_all(distr_all > q_d_all(1) & distr_all < q_d_all(2));
%     % x values for plot
%     x_val = linspace(min(d_all_show), max(d_all_show));
% %     x_val = logspace(log10(-min(d_all_show)), log10(-max(d_all_show)));
%     pdca_all = fitdist(d_all_show, 'Kernel');    
%     all_pdf = pdf(pdca_all, x_val);
%     
%           han = figure('Units', 'centimeters', ...
%                    'Position', [1, 1, 16, 20], ...
%                    'PaperSize', [16, 20]);
%     
%     % model loop
%     for m = 1:nModel
%       % values of covariance and sample set for distribution
%       distr_cov = ks_res.ReorderedTable(~isnan(ks_res.ReorderedTable(:, c+14+(r-2)*58)) & ...
%                                    ks_res.Best(:, m), c+14+(r-2)*58);
%     
%       q_d_cov = quantile(distr_cov, q_bound);  
%       % range
%       d_cov_show = distr_cov(distr_cov > q_d_cov(1) & distr_cov < q_d_cov(2));
%     
%       % fit probability distribution
%       pdca_cov = fitdist(d_cov_show, 'Kernel');    
%       cov_pdf = pdf(pdca_cov, x_val);
%     
% %       han(m) = figure('PaperSize', [14, 12]);
%       subplot(nModel/2, 2, m) 
%       area(x_val, all_pdf, 'LineWidth', 2, ...
%                            'FaceColor', 'r', ...
%                            'EdgeColor', 'r', ...
%                            'FaceAlpha', 0.2 ...
%           )
%       hold on
%       area(x_val, cov_pdf, 'LineWidth', 2, ...
%                            'FaceColor', 'b', ...
%                            'EdgeColor', 'b', ...
%                            'FaceAlpha', 0.2 ...
%           )
%       gca_act = gca;
%       axis([min(d_all_show), max(d_all_show), 0, 0.25])
%       if mod(m, 2) == 1
%         ylabel('PDF', 'Interpreter', 'latex')
%       end
% %       title([ks_res.MetafeatureNames{r, c}, ' for ', modelLabels{m}])
%       title(modelLabels{m}, 'Interpreter', 'latex')
%       legend({'all', modelLabels{m}}, 'Interpreter', 'latex')
%       hold off
%     end
% %     h_cov = histfit(distr_cov(distr_cov > q_d_cov(1) & distr_cov < q_d_cov(2)), ...
% %                     100, 'kernel');
% %     h_all = histfit(distr_all(distr_all > q_d_all(1) & distr_all < q_d_all(2)), ...
% %                     100, 'kernel');
% %     figure(4)
%     
%   end
% end
% 
% % distrFigNames = cellfun(@(x) fullfile(plotResultsFolder, ['skew_', x, '.pdf']), ...
% %                         modelLabels, 'UniformOutput', false);
% distrFigNames = {fullfile(plotResultsFolder,     'archive_skewness.pdf')};
% print2pdf(han, distrFigNames, 1)

%%

% finalize script
clear
close all
