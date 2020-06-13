%% ECML 2020 article plots
% Script for analysing GP covariance function dependence on data sampled 
% from DTS-CMA-ES run over the noiseless part of the BBOB framework.
% Script also creates graphs and tables.
% 
% Created for ECML 2020 journal track article.

%% 

% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'ecml2020_data.mat');
if isfile(tmpFName)
  load(tmpFName);
else

%%

% folders for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ecml2020abstract');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
experimentFolder = fullfile('exp', 'experiments', 'ECML2020');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);
[~, ~] = mkdir(experimentFolder);

% create tables from calculated features and save them

expfolder = fullfile('exp', 'experiments');
exp_id = 'exp_DTSmeta_02';
exp_output = fullfile(expfolder, exp_id, 'meta_output');
exp_reduced_output = fullfile(expfolder, exp_id, 'dataset', 'sampled_mfts_reduced');
[~, ~] = mkdir(exp_output);
[~, ~] = mkdir(exp_reduced_output);
exp_meta_minmax = fullfile(exp_output, 'exp_DTSmeta_02_minmax.mat');
exp_meta_inf = fullfile(exp_output, 'exp_DTSmeta_02_inf.mat');
exp_meta_output = fullfile(exp_output, 'exp_DTSmeta_02_res.mat');
exp_meta_quantile = fullfile(exp_output, 'exp_DTSmeta_02_quantile.mat');
printScriptMess = false;

% file list
fileList = searchFile('exp/experiments/exp_DTSmeta_02/dataset/sampled_metafeatures', '*.mat*');
nResFiles = numel(fileList);

%%

% Metafeature names

mfts_names = {...
     'archive_basic_dim', ...
     'archive_basic_observations', ...
     'archive_basic_lower_min', ...
     'archive_basic_lower_max', ...
     'archive_basic_upper_min', ...
     'archive_basic_upper_max', ...
     'archive_basic_objective_min', ...
     'archive_basic_objective_max', ...
     'archive_basic_blocks_min', ...
     'archive_basic_blocks_max', ...
     'archive_basic_cells_total', ...
     'archive_basic_cells_filled', ...
     'archive_basic_minimize_fun', ...
     'archive_cmaes_cma_generation', ...
     'archive_cmaes_cma_step_size', ...
     'archive_cmaes_cma_restart', ...
     'archive_cmaes_cma_mean_dist', ...
     'archive_cmaes_cma_evopath_c_norm', ...
     'archive_cmaes_cma_evopath_s_norm', ...
     'archive_cmaes_cma_lik', ...
     'archive_dispersion_ratio_mean_02', ...
     'archive_dispersion_ratio_median_02', ...
     'archive_dispersion_diff_mean_02', ...
     'archive_dispersion_diff_median_02', ...
     'archive_dispersion_ratio_mean_05', ...
     'archive_dispersion_ratio_median_05', ...
     'archive_dispersion_diff_mean_05', ...
     'archive_dispersion_diff_median_05', ...
     'archive_dispersion_ratio_mean_10', ...
     'archive_dispersion_ratio_median_10', ...
     'archive_dispersion_diff_mean_10', ...
     'archive_dispersion_diff_median_10', ...
     'archive_dispersion_ratio_mean_25', ...
     'archive_dispersion_ratio_median_25', ...
     'archive_dispersion_diff_mean_25', ...
     'archive_dispersion_diff_median_25', ...
     'archive_ela_distribution_skewness', ...
     'archive_ela_distribution_kurtosis', ...
     'archive_ela_distribution_number_of_peaks', ...
     'archive_ela_levelset_mmce_lda_10', ...
     'archive_ela_levelset_mmce_qda_10', ...
     'archive_ela_levelset_mmce_mda_10', ...
     'archive_ela_levelset_lda_qda_10', ...
     'archive_ela_levelset_lda_mda_10', ...
     'archive_ela_levelset_qda_mda_10', ...
     'archive_ela_levelset_mmce_lda_25', ...
     'archive_ela_levelset_mmce_qda_25', ...
     'archive_ela_levelset_mmce_mda_25', ...
     'archive_ela_levelset_lda_qda_25', ...
     'archive_ela_levelset_lda_mda_25', ...
     'archive_ela_levelset_qda_mda_25', ...
     'archive_ela_levelset_mmce_lda_50', ...
     'archive_ela_levelset_mmce_qda_50', ...
     'archive_ela_levelset_mmce_mda_50', ...
     'archive_ela_levelset_lda_qda_50', ...
     'archive_ela_levelset_lda_mda_50', ...
     'archive_ela_levelset_qda_mda_50', ...
     'archive_ela_metamodel_lin_simple_adj_r2', ...
     'archive_ela_metamodel_lin_simple_intercept', ...
     'archive_ela_metamodel_lin_simple_coef_min', ...
     'archive_ela_metamodel_lin_simple_coef_max', ...
     'archive_ela_metamodel_lin_simple_coef_max_by_min', ...
     'archive_ela_metamodel_lin_w_interact_adj_r2', ...
     'archive_ela_metamodel_quad_simple_adj_r2', ...
     'archive_ela_metamodel_quad_simple_cond', ...
     'archive_ela_metamodel_quad_w_interact_adj_r2', ...
     'archive_infocontent_h_max', ...
     'archive_infocontent_eps_s', ...
     'archive_infocontent_eps_max', ...
     'archive_infocontent_m0', ...
     'archive_infocontent_eps_ratio', ...
     'archive_nearest_better_nb_std_ratio', ...
     'archive_nearest_better_nb_mean_ratio', ...
     'archive_nearest_better_nb_cor', ...
     'archive_nearest_better_dist_ratio', ...
     'archive_nearest_better_nb_fitness_cor', ...
     'train_cmaes_cma_generation', ...
     'train_cmaes_cma_step_size', ...
     'train_cmaes_cma_restart', ...
     'train_cmaes_cma_mean_dist', ...
     'train_cmaes_cma_evopath_c_norm', ...
     'train_cmaes_cma_evopath_s_norm', ...
     'train_cmaes_cma_lik', ...
     'train_basic_dim', ...
     'train_basic_observations', ...
     'train_basic_lower_min', ...
     'train_basic_lower_max', ...
     'train_basic_upper_min', ...
     'train_basic_upper_max', ...
     'train_basic_objective_min', ...
     'train_basic_objective_max', ...
     'train_basic_blocks_min', ...
     'train_basic_blocks_max', ...
     'train_basic_cells_total', ...
     'train_basic_cells_filled', ...
     'train_basic_minimize_fun', ...
     'train_dispersion_ratio_mean_02', ...
     'train_dispersion_ratio_median_02', ...
     'train_dispersion_diff_mean_02', ...
     'train_dispersion_diff_median_02', ...
     'train_dispersion_ratio_mean_05', ...
     'train_dispersion_ratio_median_05', ...
     'train_dispersion_diff_mean_05', ...
     'train_dispersion_diff_median_05', ...
     'train_dispersion_ratio_mean_10', ...
     'train_dispersion_ratio_median_10', ...
     'train_dispersion_diff_mean_10', ...
     'train_dispersion_diff_median_10', ...
     'train_dispersion_ratio_mean_25', ...
     'train_dispersion_ratio_median_25', ...
     'train_dispersion_diff_mean_25', ...
     'train_dispersion_diff_median_25', ...
     'train_ela_distribution_skewness', ...
     'train_ela_distribution_kurtosis', ...
     'train_ela_distribution_number_of_peaks', ...
     'train_ela_levelset_mmce_lda_10', ...
     'train_ela_levelset_mmce_qda_10', ...
     'train_ela_levelset_mmce_mda_10', ...
     'train_ela_levelset_lda_qda_10', ...
     'train_ela_levelset_lda_mda_10', ...
     'train_ela_levelset_qda_mda_10', ...
     'train_ela_levelset_mmce_lda_25', ...
     'train_ela_levelset_mmce_qda_25', ...
     'train_ela_levelset_mmce_mda_25', ...
     'train_ela_levelset_lda_qda_25', ...
     'train_ela_levelset_lda_mda_25', ...
     'train_ela_levelset_qda_mda_25', ...
     'train_ela_levelset_mmce_lda_50', ...
     'train_ela_levelset_mmce_qda_50', ...
     'train_ela_levelset_mmce_mda_50', ...
     'train_ela_levelset_lda_qda_50', ...
     'train_ela_levelset_lda_mda_50', ...
     'train_ela_levelset_qda_mda_50', ...
     'train_ela_metamodel_lin_simple_adj_r2', ...
     'train_ela_metamodel_lin_simple_intercept', ...
     'train_ela_metamodel_lin_simple_coef_min', ...
     'train_ela_metamodel_lin_simple_coef_max', ...
     'train_ela_metamodel_lin_simple_coef_max_by_min', ...
     'train_ela_metamodel_lin_w_interact_adj_r2', ...
     'train_ela_metamodel_quad_simple_adj_r2', ...
     'train_ela_metamodel_quad_simple_cond', ...
     'train_ela_metamodel_quad_w_interact_adj_r2', ...
     'train_infocontent_h_max', ...
     'train_infocontent_eps_s', ...
     'train_infocontent_eps_max', ...
     'train_infocontent_m0', ...
     'train_infocontent_eps_ratio', ...
     'train_nearest_better_nb_std_ratio', ...
     'train_nearest_better_nb_mean_ratio', ...
     'train_nearest_better_nb_cor', ...
     'train_nearest_better_dist_ratio', ...
     'train_nearest_better_nb_fitness_cor', ...
     'traintest_cmaes_cma_generation', ...
     'traintest_cmaes_cma_step_size', ...
     'traintest_cmaes_cma_restart', ...
     'traintest_cmaes_cma_mean_dist', ...
     'traintest_cmaes_cma_evopath_c_norm', ...
     'traintest_cmaes_cma_evopath_s_norm', ...
     'traintest_cmaes_cma_lik', ...
     'traintest_basic_dim', ...
     'traintest_basic_observations', ...
     'traintest_basic_lower_min', ...
     'traintest_basic_lower_max', ...
     'traintest_basic_upper_min', ...
     'traintest_basic_upper_max', ...
     'traintest_basic_objective_min', ...
     'traintest_basic_objective_max', ...
     'traintest_basic_blocks_min', ...
     'traintest_basic_blocks_max', ...
     'traintest_basic_cells_total', ...
     'traintest_basic_cells_filled', ...
     'traintest_basic_minimize_fun', ...
     'traintest_dispersion_ratio_mean_02', ...
     'traintest_dispersion_ratio_median_02', ...
     'traintest_dispersion_diff_mean_02', ...
     'traintest_dispersion_diff_median_02', ...
     'traintest_dispersion_ratio_mean_05', ...
     'traintest_dispersion_ratio_median_05', ...
     'traintest_dispersion_diff_mean_05', ...
     'traintest_dispersion_diff_median_05', ...
     'traintest_dispersion_ratio_mean_10', ...
     'traintest_dispersion_ratio_median_10', ...
     'traintest_dispersion_diff_mean_10', ...
     'traintest_dispersion_diff_median_10', ...
     'traintest_dispersion_ratio_mean_25', ...
     'traintest_dispersion_ratio_median_25', ...
     'traintest_dispersion_diff_mean_25', ...
     'traintest_dispersion_diff_median_25', ...
     'traintest_ela_distribution_skewness', ...
     'traintest_ela_distribution_kurtosis', ...
     'traintest_ela_distribution_number_of_peaks', ...
     'traintest_ela_levelset_mmce_lda_10', ...
     'traintest_ela_levelset_mmce_qda_10', ...
     'traintest_ela_levelset_mmce_mda_10', ...
     'traintest_ela_levelset_lda_qda_10', ...
     'traintest_ela_levelset_lda_mda_10', ...
     'traintest_ela_levelset_qda_mda_10', ...
     'traintest_ela_levelset_mmce_lda_25', ...
     'traintest_ela_levelset_mmce_qda_25', ...
     'traintest_ela_levelset_mmce_mda_25', ...
     'traintest_ela_levelset_lda_qda_25', ...
     'traintest_ela_levelset_lda_mda_25', ...
     'traintest_ela_levelset_qda_mda_25', ...
     'traintest_ela_levelset_mmce_lda_50', ...
     'traintest_ela_levelset_mmce_qda_50', ...
     'traintest_ela_levelset_mmce_mda_50', ...
     'traintest_ela_levelset_lda_qda_50', ...
     'traintest_ela_levelset_lda_mda_50', ...
     'traintest_ela_levelset_qda_mda_50', ...
     'traintest_infocontent_h_max', ...
     'traintest_infocontent_eps_s', ...
     'traintest_infocontent_eps_max', ...
     'traintest_infocontent_m0', ...
     'traintest_infocontent_eps_ratio', ...
     'traintest_nearest_better_nb_std_ratio', ...
     'traintest_nearest_better_nb_mean_ratio', ...
     'traintest_nearest_better_nb_cor', ...
     'traintest_nearest_better_dist_ratio', ...
     'traintest_nearest_better_nb_fitness_cor', ...
  };

%% 

% Removing metafeatures
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
% and nbc metafeatures are useless in case of traintest
rem3_mfts = {...
     'traintest_nearest_better_nb_std_ratio', ...
     'traintest_nearest_better_nb_mean_ratio', ...
     'traintest_nearest_better_nb_cor', ...
     'traintest_nearest_better_dist_ratio', ...
     'traintest_nearest_better_nb_fitness_cor', ...
  };
rem_mfts = [rem_mfts, rem2_mfts, rem3_mfts];

fprintf('Removing user-defined metafeatures:\n')
removeId = false(219, 1);
for m = 1:numel(rem_mfts)
  mId = strcmp(mfts_names, rem_mfts{m});
  fprintf('%s\n', rem_mfts{m})
  removeId(mId) = true;
end
% mfts = mfts(:, ~mId);
mfts_names = mfts_names(~removeId);
% vars = vars(~removeId, :);
% means = means(~removeId, :);

% rename archive identical columns (see above)
mfts_names{strcmp(mfts_names, 'archive_basic_dim')} = 'dimension';
mfts_names{strcmp(mfts_names, 'archive_cmaes_cma_evopath_c_norm')} = 'cmaes_evopath_c_norm';
mfts_names{strcmp(mfts_names, 'archive_cmaes_cma_evopath_s_norm')} = 'cmaes_evopath_s_norm';
mfts_names{strcmp(mfts_names, 'archive_cmaes_cma_generation')}     = 'cmaes_generation';
mfts_names{strcmp(mfts_names, 'archive_cmaes_cma_restart')}        = 'cmaes_restart';
mfts_names{strcmp(mfts_names, 'archive_cmaes_cma_step_size')}      = 'cmaes_step_size';
% rename observations columns
mfts_names{strcmp(mfts_names, 'archive_basic_observations')}   = 'archive_observations';
mfts_names{strcmp(mfts_names, 'train_basic_observations')}     = 'train_observations';
mfts_names{strcmp(mfts_names, 'traintest_basic_observations')} = 'traintest_observations';

% create mfts table, where only user defined metafeatures are not present
% full_mfts = mfts;

% metafeature names for printing
mfts_names_prt = strrep(mfts_names, '_', '\_');

clear rem_mfts rem2_mfts rem3_mfts m ms mId mss

%% 

% reduce results to feature values, names, and basic identifiers

% list of ids to reduce (exclude n.7)
redIds = [1:6, 8:9];

for fl = 1:nResFiles
  [~, actFile] = fileparts(fileList{fl});
  reducedFilename = fullfile(exp_reduced_output, [actFile, '_red.mat']);
  id_str = regexp(actFile, '_id(\d)_', 'tokens');
  id = str2double(id_str{1});
  % create reduced results
  if ~isfile(reducedFilename) && ismember(id, redIds)
    if printScriptMess
      fprintf('Saving %s\n', reducedFilename)
    end
    S = load(fileList{fl});
    mfts_values = S.res.values(~removeId, :, :);
    dim = S.dim;
    fun = S.fun;
    inst = S.inst;
    save(reducedFilename, 'dim', 'fun', 'inst', 'id', 'mfts_names', 'mfts_values')
  end
end

clear actFile fl dim fun inst id id_str mfts_values S reducedFilename

%% 

% feature distribution

% create list of files with reduced data
redFileList = searchFile(exp_reduced_output, '*.mat*');

q_levels = 0:0.01:1;

if isfile(exp_meta_quantile)
  S = load(exp_meta_quantile);
  q_ft = S.q_ft;
else
  for ft = 1:numel(mfts_names)
    fprintf('%s\n', mfts_names{ft})
    
    actual_ft = [];
    for fl = 1: numel(redFileList)
      S = load(redFileList{fl}, 'mfts_values');
      actual_ft = [actual_ft, shiftdim(S.mfts_values(ft, :, :), 1)];
    end
    q_ft(ft, :) = quantile(actual_ft(:), q_levels);
  end
  save(exp_meta_quantile, 'q_ft')
end

clear actual_ft fl ft q_levels S

%% 

% feature mins and maxs

% create list of files with reduced data
redFileList = searchFile(exp_reduced_output, '*.mat*');

if isfile(exp_meta_minmax)
  MM = load(exp_meta_minmax);
  ft_max = MM.ft_max;
  ft_min = MM.ft_min;
  ft_min_ninf = MM.ft_min_ninf;
  ft_max_ninf = MM.ft_max_ninf;
else
  if printScriptMess
    fprintf('Calculating feature minimums and maximums') 
  end
  S = load(exp_meta_quantile, 'q_ft');
  ft_min = S.q_ft(:, 1);
  ft_max = S.q_ft(:, end);
  
  % find inf valued min and max
  ft_inf = isinf(ft_min) | isinf(ft_max);
  
  % init values
  ft_min_ninf = ft_min;
  ft_max_ninf = ft_max;
  ft_min_ninf(ft_inf) = NaN(sum(ft_inf), 1);
  ft_max_ninf(ft_inf) = NaN(sum(ft_inf), 1);
  
  % compare with the rest of results
  for fl = 1:numel(redFileList)
    S = load(redFileList{fl}, 'mfts_values');
    % find non-inf min and max values
    for ft = find(ft_inf)'
      actMfts = S.mfts_values(ft, :, :);
      act_ft_min = min(actMfts(~isinf(actMfts)));
      act_ft_max = max(actMfts(~isinf(actMfts)));
      if isempty(act_ft_min)
        S_ft_min_ninf(ft, 1) = NaN;
      else
        S_ft_min_ninf(ft, 1) = act_ft_min;
      end
      if isempty(act_ft_max)
        S_ft_max_ninf(ft, 1) = NaN;
      else
        S_ft_max_ninf(ft, 1) = act_ft_max;
      end
    end
    % compare to overall 2nd min and max
    ft_max_ninf = max(ft_max_ninf, S_ft_max_ninf);
    ft_min_ninf = min(ft_min_ninf, S_ft_min_ninf);
  end
  
  % save inf and non-inf min and max
  save(exp_meta_minmax, 'ft_min', 'ft_max', 'ft_min_ninf', 'ft_max_ninf')
end

clear actMfts act_ft_max act_ft_min fl ft MM S S_ft_max_ninf S_ft_min_ninf

%%

% feature infs

% create list of files with reduced data
redFileList = searchFile(exp_reduced_output, '*.mat*');

if isfile(exp_meta_inf)
  MM = load(exp_meta_inf);
  inf_plus = MM.inf_plus;
  inf_min = MM.inf_min;
else
  if printScriptMess
    fprintf('Calculating feature infs') 
  end
  S = load(redFileList{1}, 'mfts_values');
  % calculate number of infs
  inf_plus = sum(sum(isinf(S.mfts_values) & S.mfts_values > 0, 3), 2);
  inf_min = sum(sum(isinf(S.mfts_values) & S.mfts_values < 0, 3), 2);
  
  % compare with the rest of results
  for fl = 2:numel(redFileList)
    S = load(redFileList{fl}, 'mfts_values');
    % calculate number of infs
    inf_plus = inf_plus + sum(sum(isinf(S.mfts_values) & S.mfts_values > 0, 3), 2);
    inf_min = inf_min + sum(sum(isinf(S.mfts_values) & S.mfts_values < 0, 3), 2);
  end
  
  % save min and max
  save(exp_meta_inf, 'inf_plus', 'inf_min')
end

clear fl MM S

%%

% load and normalize results

% create list of files with reduced data
redFileList = searchFile(exp_reduced_output, '*.mat*');

if isfile(exp_meta_minmax)
  load(exp_meta_minmax, 'ft_min', 'ft_max', 'ft_min_ninf', 'ft_max_ninf')
else
  error('There is %s missing.', exp_meta_minmax)
end

if ~isfile(exp_meta_output)
  if printScriptMess
    fprintf('Loading metalearing results\n')
  end
  
  % init  
  observations = ([]); % expected range [0, 5000]
  generations = ([]);  % expected range [0, 5000]
  dimensions = ([]);    % expected range [2, 20]
  vars = [];
  means = [];
  medians = [];
  nans = ([]);     % expected range [0, 100]
  inf_plus = ([]); % expected range [0, 100]
  inf_mins = ([]); % expected range [0, 100]

  for fl = 1:numel(redFileList)
    S = load(redFileList{fl}, 'mfts_values');

    % normalize to [0,1]
    for ft = 1:size(S.mfts_values, 1)
      act_mfts_values = S.mfts_values(ft, :, :);
      
      % linear transformation to [0, 1]
      norm_val(ft, :, :) = (act_mfts_values - ft_min_ninf(ft))/(ft_max_ninf(ft)-ft_min_ninf(ft));
    end

    % collect generations, numbers of points, and dimensions
    observations = [observations, [S.mfts_values(2, :, 1); S.mfts_values(67, :, 1); S.mfts_values(125, :, 1)]];
    generations  = [generations,  S.mfts_values(3, :, 1)];
    dimensions   = [dimensions,   S.mfts_values(1, :, 1)];
    
    % collect nans and infs
    nans = [nans, sum(isnan(norm_val), 3)];
    inf_plus = [inf_plus, sum(isinf(norm_val) & norm_val > 0, 3)];
    inf_mins = [inf_mins, sum(isinf(norm_val) & norm_val < 0, 3)];
    
    % calculate statistics from non-inf and non-nan values
    norm_val(isinf(norm_val)) = NaN;
    vars  = [vars,  nanvar(norm_val, 0, 3)];
    means = [means, nanmean(norm_val, 3)];
    medians = [medians, nanmedian(norm_val, 3)];
    
    % save results after each 100 files
    if mod(fl, 100) == 0
      save(exp_meta_output, 'observations', 'generations', 'dimensions', 'nans', 'inf_plus', 'inf_mins', 'vars', 'means', 'medians')
    end
  end
  
  % save overall results
  save(exp_meta_output, 'observations', 'generations', 'dimensions', 'nans', 'inf_plus', 'inf_mins', 'vars', 'means', 'medians')
%   save(exp_meta_output, 'nans', 'inf_plus', 'inf_mins', 'vars', 'means', 'medians')
end

clear actualFolder articleFolder fl ft i infBound logscale norm_val S act_ft_min act_ft_max act_mfts_values

%%

if isfile(exp_meta_output)
  load(exp_meta_output)
else
  error('File %s is missing!', exp_meta_output)
end

% overall results
if printScriptMess
  fprintf('Overall results\n')
end
overall_var  = nanmedian(vars, 2);
overall_mean = nanmedian(means, 2);

% results according to dimension
if printScriptMess
  fprintf('Dimension results\n')
end
dims = [2, 3, 5, 10, 20];
for d = 1:numel(dims)
  dim_var(:, d)  = nanmedian(vars(:, dims(d) == dimensions), 2);
  dim_mean(:, d) = nanmedian(means(:, dims(d) == dimensions), 2);
end

% number of xaxis quantiles
nQuant = 1000;
nMfts = numel(mfts_names);
% plot data quantiles
disp_quant = [0.05, 0.5, 0.95];
nDispQuant = numel(disp_quant);

% init variables
obs_vars = NaN(nMfts, nQuant, nDispQuant);
obs_means = NaN(nMfts, nQuant, nDispQuant);
obs_nans = NaN(nMfts, nQuant, nDispQuant);
obs_inf_plus = NaN(nMfts, nQuant, nDispQuant);
obs_inf_mins = NaN(nMfts, nQuant, nDispQuant);
obs_dims = false(size(observations, 1), nQuant, numel(dims));

dens_vars = NaN(nMfts, nQuant, nDispQuant);
dens_means = NaN(nMfts, nQuant, nDispQuant);
dens_nans = NaN(nMfts, nQuant, nDispQuant);
dens_inf_plus = NaN(nMfts, nQuant, nDispQuant);
dens_inf_mins = NaN(nMfts, nQuant, nDispQuant);
dens_dims = false(size(observations, 1), nQuant, numel(dims));

if printScriptMess
  fprintf('Results according to the number of observations\n')
end
sampleSets = {'archive_', 'train_', 'traintest_'};
% calculate variance and mean medians according to number of observations 
% for different sample sets
noSampleSetId = ~contains(mfts_names, sampleSets);
nGenFeat = sum(noSampleSetId); % number of sample independent features
un_observations = cell(1, 3);
for o = 1:size(observations, 1)
  un_observations{o} = unique(observations(o, :));
  sampleSetId = contains(mfts_names, sampleSets{o});
  % calculate densities
  densities(o, :) = nthroot(observations(o, :), dimensions);
  % get ids of sorted values
  [~,  obs_I] = sort(observations(o, :));
  [~, dens_I] = sort(   densities(o, :));
  % exclude zero observations
  obs_I(~isnatural(observations(o, obs_I))) = [];
  dens_I(~(densities(o, dens_I) > 0)) = [];
  % create permilles (1000-quantile) data groups
  obs_gId  = ceil(nQuant*(1:numel(obs_I))  / numel(obs_I));
  dens_gId = ceil(nQuant*(1:numel(dens_I)) / numel(dens_I));
  
  for q = 1:nQuant
    % observations
    act_dataId = obs_I(obs_gId == q);
    obs_vars(sampleSetId, q, :)      = shiftdim(quantile(     vars(sampleSetId, act_dataId), disp_quant, 2), -1);
    obs_means(sampleSetId, q, :)     = shiftdim(quantile(    means(sampleSetId, act_dataId), disp_quant, 2), -1);
    obs_nans(sampleSetId, q, :)      = shiftdim(quantile(     nans(sampleSetId, act_dataId), disp_quant, 2), -1);
    obs_inf_plus(sampleSetId, q, :)  = shiftdim(quantile( inf_plus(sampleSetId, act_dataId), disp_quant, 2), -1);
    obs_inf_mins(sampleSetId, q, :)  = shiftdim(quantile( inf_mins(sampleSetId, act_dataId), disp_quant, 2), -1);
    % dimensions present in particular quantile
    obs_dims(o, q, :) = ismember(dims, dimensions(act_dataId));
    % observations in particular quantile
    obs_obs(o, q, :) = median(observations(o, act_dataId));
    % calculate only means for sample independent features
%     general_obs_means((1:nGenFeat) + (o-1)*nGenFeat, uo) = ...
%       median(means(noSampleSetId, un_observations{o}(uo) == observations(o, :)), 2);

    % densities
    act_dataId = dens_I(dens_gId == q);
    dens_vars(sampleSetId, q, :)      = shiftdim(quantile(     vars(sampleSetId, act_dataId), disp_quant, 2), -1);
    dens_means(sampleSetId, q, :)     = shiftdim(quantile(    means(sampleSetId, act_dataId), disp_quant, 2), -1);
    dens_nans(sampleSetId, q, :)      = shiftdim(quantile(     nans(sampleSetId, act_dataId), disp_quant, 2), -1);
    dens_inf_plus(sampleSetId, q, :)  = shiftdim(quantile( inf_plus(sampleSetId, act_dataId), disp_quant, 2), -1);
    dens_inf_mins(sampleSetId, q, :)  = shiftdim(quantile( inf_mins(sampleSetId, act_dataId), disp_quant, 2), -1);
    % dimensions present in particular quantile
    dens_dims(o, q, :) = ismember(dims, dimensions(act_dataId));
    % densities in particular quantile
    dens_dens(o, q, :) = median(densities(o, act_dataId));
  end
end

clear d o q
% clear large variables saved in exp_meta_output
clear means medians vars nans inf_plus inf_mins

%%
if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

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

% archive_dispersion_diff_median_02 -> ft  13
% traintest_infocontent_eps_s       -> ft 164

h = [];
figNames = {};

close all

for ft = [13, 164] %1:numel(mfts_names)
  %%
  
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
%     mftPropPlot(shiftdim(obs_means(ft, :, :), 1), med_obs_nans, ...
%                 med_obs_inf_plus, med_obs_inf_mins, ...
%                 'Dimensions', shiftdim(obs_dims(sampleSetId(ft), :, :), 1), ...
%                 'DimTypes', dims, ...
%                 'MftName', mfts_names_prt{ft}, ...
%                 'ValueName', {'0.05 mean', '0.5 mean', '0.95 mean'}, ...
%                 'ValueStyle', {'-.', '-', '-.'}, ...
%                 'XLabelName', sprintf('%d-quantile of observations', nQuant));
              
    h(end+1) = mftPropPlot(shiftdim(obs_vars(ft, :, :), 1), med_obs_nans, ...
                med_obs_inf_plus, med_obs_inf_mins, ...
                'Dimensions', [], ... shiftdim(obs_dims(sampleSetId(ft), :, :), 1), ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 var', '0.5 var', '0.95 var'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', obs_obs(sampleSetId(ft), :), ...
                'XLabelName', 'number of observations');
    figNames{end+1} = fullfile(plotResultsFolder, [mfts_names{ft}, '_var_obs.pdf']);
              
    % densities
%     mftPropPlot(shiftdim(dens_means(ft, :, :), 1), med_dens_nans, ...
%                 med_dens_inf_plus, med_dens_inf_mins, ...
%                 'Dimensions', shiftdim(dens_dims(sampleSetId(ft), :, :), 1), ...
%                 'DimTypes', dims, ...
%                 'MftName', mfts_names_prt{ft}, ...
%                 'ValueName', {'0.05 mean', '0.5 mean', '0.95 mean'}, ...
%                 'ValueStyle', {'-.', '-', '-.'}, ...
%                 'XLabelName', sprintf('%d-quantile of densities', nQuant));
              
    h(end+1) = mftPropPlot(shiftdim(dens_vars(ft, :, :), 1), med_dens_nans, ...
                med_dens_inf_plus, med_dens_inf_mins, ...
                'Dimensions', [], ...
                'MftName', mfts_names_prt{ft}, ...
                'ValueName', {'0.05 var', '0.5 var', '0.95 var'}, ...
                'ValueStyle', {'-.', '-', '-.'}, ...
                'XValues', dens_dens(sampleSetId(ft), :), ...
                'XValuesDec', 2, ...
                'XLabelName', 'density of observations');
    figNames{end+1} = fullfile(plotResultsFolder, [mfts_names{ft}, '_var_dens.pdf']);

  end
end

print2pdf(h, figNames, 1);

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
