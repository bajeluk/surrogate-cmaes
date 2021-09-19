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
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
experimentFolder = fullfile('exp', 'experiments', 'exp_smsp');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);
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
mfts_names{strcmp(mfts_names, [identicalLabel, 'basic_dim'])} = 'dimension';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_evopath_c_norm'])} = 'cmaes_evopath_c_norm';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_evopath_s_norm'])} = 'cmaes_evopath_s_norm';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_generation'])}     = 'cmaes_generation';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_restart'])}        = 'cmaes_restart';
mfts_names{strcmp(mfts_names, [identicalLabel, 'cmaes_cma_step_size'])}      = 'cmaes_step_size';
% rename observations columns
for m = 1:numel(mfts_names)
  if contains(mfts_names{m}, 'observations')
    mfts_name_parts = strsplit(mfts_names{m}, '_');
    mfts_names{m} = strjoin({mfts_name_parts{1}, 'observations'}, '_');
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
      observations(:, colId) = S.mfts_values(contains(mfts_names_red, 'observations'), :, 1);
      if tss == fullId
        generations(:, colId) = S.mfts_values(contains(mfts_names_red, 'generation'),   :, 1);
        dimensions(:, colId)  = S.mfts_values(contains(mfts_names_red, 'dimension'),    :, 1);
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

  % load results
  if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'nans', 'inf_plus', 'inf_mins', 'mfts_names_red');
  else
    error('There is %s missing.', exp_meta_output{tss})
  end
  mfts_names_red = S.mfts_names_red;

  nan_perc = sum(S.nans, 2)/size(S.nans, 2);
  inf_perc = sum(S.inf_mins+S.inf_plus, 2)/size(S.nans, 2);

  % save stats
  save(exp_meta_stats{tss}, 'inf_perc', 'nan_perc', 'mfts_names_red')
  
  % mark features with high NaN percentage
  mftsLowNaN = nan_perc < 25;
  save(exp_meta_output{tss}, 'mftsLowNaN', '-append')
end

clear inf_perc nan_perc S tss

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
    load(exp_meta_output{tss})
  else
    error('File %s is missing!', exp_meta_output{tss})
  end
  
  % use only features with low number of NaNs
  mfts_names_red = mfts_names_red(mftsLowNaN);
  medians = medians(mftsLowNaN, :);

  nMfts = numel(mfts_names_red);
  % medtest_meds_p = zeros(nMfts, nchoosek(numel(dims), 2));
  friedman_meds_p = NaN(nMfts, 1);
  stats_meds = cell(nMfts, 1);
  multcomp_meds = cell(nMfts, 1);
  friedman_meds_p_bh = NaN(nMfts, nchoosek(numel(dims), 2));
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
    dimSizes = arrayfun(@(x) sum(x == dimensions), dims);
    minDimSize = min(dimSizes);
    if any(dimSizes ~= minDimSize)
      warning('The size of results in different dimensions are not the same for %s. Reducing to the lowest size.', tssList{tss})
    end
    % cat medians for individual dimensions (number of data for each
    % dimension must be identical)
    for d = 1:numel(dims)
      mat_meds_act = act_meds(dims(d) == dimensions)';
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
          meds_pair = [act_meds(dims(d) == dimensions); ...
                       act_meds(dims(d2) == dimensions)];
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

  if isfile(exp_meta_output{tss})
    S = load(exp_meta_output{tss}, 'observations', ...
             'nans', 'mftsLowNaN', ...
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
  mfts_names_red = S.mfts_names_red(S.mftsLowNaN);
  S.nans = S.nans(S.mftsLowNaN, :);
  
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
             'quantiles_90', 'quantiles_99', 'quantile_to_calc', ...
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
  mfts_names_red = S.mfts_names_red(S.mftsLowNaN);
  S.quantiles_90 = S.quantiles_90(S.mftsLowNaN, :, :);
  S.quantiles_99 = S.quantiles_99(S.mftsLowNaN, :, :);
  
  % init
  nMfts = numel(mfts_names_red);
  quan_90_diff_break = zeros(nMfts, 1);
  quan_90_diff_meds = cell(nMfts, 1);
  quan_99_diff_break = zeros(nMfts, 1);
  quan_99_diff_meds = cell(nMfts, 1);

  % feature loop
  for ft = 1:nMfts

    %%
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
      quan_90_diff_meds{ft}(i) = median(...
        S.quantiles_90(ft, uniObsId, 2) - S.quantiles_90(ft, uniObsId, 1));
      quan_99_diff_meds{ft}(i) = median(...
        S.quantiles_99(ft, uniObsId, 2) - S.quantiles_99(ft, uniObsId, 1));
    end

    % plot resulting data
    figure()
    plot(obsUni{obsId}, quan_90_diff_meds{ft})
    hold on
    plot(obsUni{obsId}, quan_99_diff_meds{ft})
    line([0, obsUni{obsId}(end)], [reliability_threshold, reliability_threshold], 'Color', 'r')
    title(['TSS ', num2str(tss), ': ', strrep(mfts_names_red{ft}, '_', '\_')])
    legend('Quantile diff (90%)', 'Quantile diff (99%)', '0.05 reliability threshold')
    hold off
  end

  % save results
  save(exp_smsp_quant_anal{tss}, 'nan_break', 'nan_meds', 'nan_threshold', ...
                               'mfts_names_red', 'obsUni')
end

clear('ft', 'i', 'mfts_names_red', 'nan_break', 'nan_breakId', ...
      'nan_meds', 'nMfts', 'obsId', 'obsUni', 'tss')

close

%% Feature correlation analysis
% Schweizer-Wolff correlation between each pair of available features used
% to calculation of distances between features

nPerms = 5;
threshold = 0.25;

for tss = 1:nTSS
  fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Performing feature correlation analysis on TSS %s\n'], ...
          fix(clock), tssList{tss})
  % for each TSS load data for TSS full
  if isfile(exp_meta_output{1})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS full output\n'], ...
          fix(clock), tssList{tss})
    S_full = load(exp_meta_output{1}, ...
                  'medians', 'inf_mins', 'inf_plus', 'mfts_names_red', ...
                  'mftsLowNaN');
    % use only features with low number of NaNs
    S_full.mfts_names_red = S_full.mfts_names_red(S_full.mftsLowNaN);
    S_full.medians = S_full.medians(S_full.mftsLowNaN, :);
    S_full.inf_mins = S_full.inf_mins(S_full.mftsLowNaN, :);
    S_full.inf_plus = S_full.inf_plus(S_full.mftsLowNaN, :);
  else
    error('File %s is missing!', exp_meta_output{1})
  end
  % load minmax TSS full
  if isfile(exp_meta_minmax{1})
    fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS full minmax\n'], ...
          fix(clock), tssList{tss})
    S_mm_full = load(exp_meta_minmax{1}, 'ft_min_ninf', 'ft_max_ninf');
    % use only features with low number of NaNs
    S_mm_full.ft_min_ninf = S_mm_full.ft_min_ninf(S_full.mftsLowNaN);
    S_mm_full.ft_max_ninf = S_mm_full.ft_max_ninf(S_full.mftsLowNaN);
  else
    error('File %s is missing!', exp_meta_minmax{1})
  end

  if tss > 1
    % load output
    if isfile(exp_meta_output{tss})
      fprintf(['[ %d-%.2d-%.2d %.2d:%.2d:%.2d ] ', ...
           'Correlation analysis on TSS %s - loading TSS %s output\n'], ...
          fix(clock), tssList{tss}, tssList{tss})
      S_tss = load(exp_meta_output{tss}, ...
                  'medians', 'inf_mins', 'inf_plus', 'mfts_names_red', ...
                  'mftsLowNaN');
      % use only features with low number of NaNs
      S_tss.mfts_names_red = S_tss.mfts_names_red(S_tss.mftsLowNaN);
      S_tss.medians = S_tss.medians(S_tss.mftsLowNaN, :);
      S_tss.inf_mins = S_tss.inf_mins(S_tss.mftsLowNaN, :);
      S_tss.inf_plus = S_tss.inf_plus(S_tss.mftsLowNaN, :);
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
      S_mm_tss.ft_min_ninf = S_mm_tss.ft_min_ninf(S_tss.mftsLowNaN);
      S_mm_tss.ft_max_ninf = S_mm_tss.ft_max_ninf(S_tss.mftsLowNaN);
    else
      error('File %s is missing!', exp_meta_minmax{tss})
    end
  else
    S_tss.medians = [];
    S_tss.inf_mins = [];
    S_tss.inf_plus = [];
    S_tss.mfts_names_red = {};
    S_tss.generations = [];
    S_tss.mftsLowNaN = [];
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
  mftsLowNaN = [S_full.mftsLowNaN; S_tss.mftsLowNaN];
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
                           'mftsLowNaN', ...
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
  mftsIds_red = mftsIds_red(logical(S.mftsLowNaN));
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
      nan_break = [S_nan_full.nan_break(:); S_nan_tss.nan_break(:)];
    else
      error('File %s is missing!', exp_smsp_nan_anal{tss})
    end
  else
    nan_break = S_nan_full.nan_break(:);
  end

  % TODO: dimension test was probably performed without TSS full and actual
  % TSS unification => check and perhaps unify here
  % TSS full dimension test result
  if isfile(exp_smsp_dimension_test{fullId})
    S_dim_full = load(exp_smsp_dimension_test{fullId});
  else
    error('File %s is missing!', exp_smsp_dimension_test{fullId})
  end
  % add actual TSS results
  if tss ~= fullId
    if isfile(exp_smsp_dimension_test{tss})
      S_dim_tss = load(exp_smsp_dimension_test{tss});
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
  nMfts = numel(S_corr.mfts_names_red);
  
  % dimension test results
  dimCombs = nchoosek(dims, 2);
  dimCombsStr = arrayfun(@(x) ...
                  sprintf('(%d, %d)', dimCombs(x, 1), dimCombs(x, 2)), ...
                  1:size(dimCombs, 1), 'uni', false);
  % TSS full combinations
  dimCombsAccepted = arrayfun(@(x) ...
    strjoin( ...
      dimCombsStr( ...
        S_dim_full.wilcoxon_uniq_combs( ...
          S_dim_full.wilcoxon_uniq_combs_id(x), ...
          :) ...
        ), ...
      ',' ...
    ), ...
    1:size(S_dim_full.wilcoxon_uniq_combs_id, 1), ...
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
        ',' ...
      ), ...
      1:size(S_dim_tss.wilcoxon_uniq_combs_id, 1), ...
      'UniformOutput', false)' ...
    ];

  % sort metafeatures to clusters
  [~, sortClId] = sort(S_corr.corrClusterId_nanpair);

  % is medoid id
  isMed = ismember((1:nMfts)', S_corr.corrMedoidId_nanpair);

  % table mfts notation
  mftsSplit = cellfun(@(x) strsplit(x, '_'), S_corr.mfts_names_red, 'Uni', false);
  setColumn   = cell(nMfts, 1);
  transColumn = cell(nMfts, 1);
  classColumn = cell(nMfts, 1);
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
           classColumn{m} = '\featCMA';
         case 'dispersion'
           classColumn{m} = '\featDisp';
         case 'ela'
           switch mftsSplit{m}{2}
             case 'distribution'
               classColumn{m} = '\featyDis';
             case 'levelset'
               classColumn{m} = '\featLevel';
             case 'metamodel'
               classColumn{m} = '\featMM';
           end
           % remove ela set notation
           mftsSplit{m}(1) = [];
         case 'infocontent'
           classColumn{m} = '\featInfo';
         case 'nearest'
           classColumn{m} = '\featNBC';
           % remove 'nearest' part of set notation
           mftsSplit{m}(1) = [];
       end
      % remove set notation
      mftsSplit{m}(1) = [];
    end
    % unite the rest
    mftsSplit{m} = strjoin(mftsSplit{m}, '\\_');
  end
  tableMftsNotation = cellfun(...
    @(x, y, z, w) sprintf('$\\tableFeat{\\texttt{%s}}{%s}{%s}{%s}$', x, y, z, w), ...
    mftsSplit, setColumn, transColumn, classColumn, 'Uni', false);

  % sort table variables
  nan_break = nan_break(sortClId);
  dimCombsAccSorted = dimCombsAccepted(sortClId);
  tableMftsNotationSorted = tableMftsNotation(sortClId);
  corrClusterId_nanpairSorted = S_corr.corrClusterId_nanpair(sortClId);
  isMedSorted = isMed(sortClId);
  % page cycle
  pages = ceil(nMfts/maxFtsPerPage);
  for p = 1:pages
    rowIds = (p-1)*maxFtsPerPage+1 : min(p*maxFtsPerPage, nMfts);
    % table with correlation clustering ids
    corrClusterTable = table(...
                             nan_break(rowIds), ...
                             dimCombsAccSorted(rowIds), ...
                             'RowNames', tableMftsNotationSorted(rowIds));

    % print table to tex
    lt = LatexTable(corrClusterTable);
    lt.opts.tableColumnAlignment = num2cell('lrl');
    lt.setHeaderRow({'', '$N_\nanout$', '$(\dm_i, \dm_j)$'});
    lt.setColumnFormat({'%d', []});
    [~, lt.opts.midLines] = unique(corrClusterId_nanpairSorted(rowIds));
    lt.opts.tableCaption = sprintf([...
      'TSS %s feature properties (%d/%d). ', ...
      'Features are grouped to %d clusters according to k-medoid clustering ', ...
      'using Schweizer-Wolf correlation distance. ', ...
      'Medoid representatives are marked as gray lines in clusters divided by horizontal lines. ', ...
      'If no gray line is present, the cluster is divided between two tables. ', ...
      '$N_\\nanout$ denotes the lowest measured number of points from which 0.01 of feature calculations resulted in $\\nanout$ at maximum. ', ...
      'The $(\\dm_i, \\dm_j)$ column shows the pairs of feature dimensions not rejecting the independence of median feature values, ', ...
      '\\ie where the two-sided Wilcoxon signed rank test on median values with the Bonferroni-Holm correction ', ...
      'at the family-wise level 0.05 was not rejected.'...
      ], tssList{tss}, p, pages, S_corr.k);
    lt.opts.tableLabel = 'featProp';
    lt.opts.booktabs = 1;
    % add gray background to medoid features
    lt.colorizeRowsInGray(isMedSorted(rowIds));
    lt.toFile(sprintf('%s_%d.tex', exp_smsp_corr_dim_table{tss}, p));
  end
  
  % save mfts notation
  save(exp_smsp_corr_cluster{tss}, 'tableMftsNotation', '-append')

  % show tables with cluster differences
  % diffIds = false(size(corrClusterId_all));
  % for c = 1:k
  %   if (range(corrClusterId_all(corrClusterId_pair==c)) > 0)
  %     diffIds = diffIds | corrClusterId_pair==c;
  %     corrClusterTable(corrClusterId_pair==c, :)
  %   end
  % end
end

clear('c', 'm', 'nMfts', 'S_corr', 'S_dim_full', 'S_dim_tss')

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
