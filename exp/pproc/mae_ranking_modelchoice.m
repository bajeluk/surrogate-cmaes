%% load GP model results
load('exp/experiments/exp_GPtest_01/modelStatistics.mat');

% NOISY modelStatistics
% load('exp/experiments/exp_GPtest_02_noisy/modelStatistics_noisy.mat');

%% Core settings
% dimensions = [2];

for dimensions = [2, 5, 10]

snapshots = [3, 9];
multiFieldNames = { 'covFcn', 'trainsetType', 'trainRange', 'trainsetSizeMax', 'meanFcn' };

%% Model choice settings

disp('##################################');
disp('for noisy, use:');
disp('load(''exp/experiments/exp_GPtest_02_noisy/modelStatistics_noisy.mat'');');
disp('maeRankingOpts.maxRank = {10, 20};');
disp('##################################');

% maximal allowed rank for choosing the model for the function/snapshot
maeRankingOpts.maxRank = {25, 35};
% index of the statistic by which the ranking is done, 1 = mean, 2 = 75%-quantile
maeRankingOpts.colStat = {2};
% minimal allowed training success rate for choosing the model for the function/snapshot
maeRankingOpts.minTrainedPerc = {0.70, 0.85};
% index of the statistic by which the error-ranking is done, 1 = mean, 2 = 75%-quantile
maeRankingOpts.colSetCover = {1, 2};
% transformation of the number of covered functions/snapshots
% ('fsCovered' is a column vector of these number for each settings)
% and ranks of respective models for covering functions/snapshots
% ('mRanks' is a matrix of size 'nSettings x sum(~isCovered)')
maeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
    @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };

maeRankingOpts.includeARD = { false };
maeRankingOpts.include5dim = { true };
maeRankingOpts.includeMeanLinear = { false };

% get the numbers of factors
optsFields = fieldnames(maeRankingOpts);
nFactors = zeros(size(optsFields));
for f = 1:length(optsFields)
  thisField = optsFields{f};
  nFactors(f) = length(maeRankingOpts.(thisField));
end

% create full factorial design (of the indices)
ffDesign = fullfact(nFactors);
n_fullfact = size(ffDesign, 1);

% go through all the fullfact. design
for ffi = 1:n_fullfact
  % get the actual opts
  opts = struct();
  for f = 1:length(optsFields)
    thisField = optsFields{f};
    opts.(thisField) = maeRankingOpts.(thisField){ffDesign(ffi, f)};
  end
  
  % perform the "test"
  run('mae_ranking');
  bestSettingsTables{ffi} = bestSettingsTable{1};
  bestSettingsCells{ffi} = bestSettingsCell{1};
  bestSettingsStructs{ffi} = bestSettingsStruct{1};
  bestSettingsOpts{ffi} = opts;
end

%% make a table with results
maeRankingResults_header = ['id', optsFields', ...
  { 'nSettings', 'avg_nCovered', ...
    'avg_MAE', 'avg_covrd_MAE', 'max_covrd_MAE', 'avg_rank', 'avg_covrd_rank', ...
    'avg_n_covrd_funs', 'n_ARD', 'n_Linear' }];
maeRankingResults = table((1:n_fullfact)');
for f = 1:length(optsFields)
  maeRankingResults(:,1+f) = cell2table( ...
       cellfun(@(x) x.(optsFields{f}), bestSettingsOpts, 'UniformOutput', false)');
end

for ffi = 1:n_fullfact
  % nSettings
  f = length(optsFields) + 2;
  maeRankingResults{ffi, f} = size(bestSettingsTables{ffi}, 1);
  % avg_nCovered
  f = f + 1;
  maeRankingResults{ffi, f} = nanmean(bestSettingsTables{ffi}.No_covered);
  % avg_MAE
  f = f + 1;
  maeRankingResults{ffi, f} = nanmean(bestSettingsTables{ffi}.avg_MAE);
  % avg_covrd_MAE
  f = f + 1;
  maeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.covrd_MAE);
  % max_covrd_MAE
  f = f + 1;
  maeRankingResults{ffi, f} = max(bestSettingsTables{ffi}.covrd_MAE);
  % avg_rank
  f = f + 1;
  maeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.avg_rank);
  % avg_covrd_rank
  f = f + 1;
  maeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.covrd_rank);
  % avg_n_covrd_funs
  f = f + 1;
  maeRankingResults{ffi, f} = mean( cellfun( ...
      @(x) length(strsplit(x, ' ')), bestSettingsTables{ffi}.covrd_funs) );
  % n_ARD
  f = f + 1;
  maeRankingResults{ffi, f} = sum( cellfun( @(x) ~isempty(x), ...
      regexpi(bestSettingsTables{ffi}.covFcn, 'ard') ));
  % n_Linear
  f = f + 1;
  maeRankingResults{ffi, f} = sum( cellfun( @(x) ~isempty(x), ...
      regexpi(bestSettingsTables{ffi}.meanFcn, 'Linear') ));
end
maeRankingResults.Properties.VariableNames = maeRankingResults_header;

%% output results
fprintf('\n==== %d D ==== \n\n', dimensions);
fprintf(['\nBest sets of settings according to the average MAE\n'...
    'of covered fuctions/snapshot:\n\n']);
tmp = sortrows(maeRankingResults, 'avg_covrd_MAE');
disp(tmp(1:3,:));
% disp(sortrows(maeRankingResults, 'max_covrd_MAE'));

% print the 3 best settings options structs
N_PRINT = 1;
[~, idx] = sort(maeRankingResults.avg_covrd_MAE);
% [~, idx] = sort(maeRankingResults.avg_covrd_MAE);
for i = 1:N_PRINT
  disp(bestSettingsTables{idx(i)});
  fprintf('\n== Best set of settings #%d (id = %d)== \n\n', i, idx(i));
  printStructure(bestSettingsStructs{idx(i)});
end

end  % for dimensions = [2, 5, 10]

%% find aggMAE_table rows which correspond to the historically most used settings
%{
origSettingsRows = strcmp(aggMAE_table.covFcn, '{@covMaterniso, 5}') & ...
    strcmp(aggMAE_table.trainsetType, 'clustering') & ...
    aggMAE_table.trainRange == 0.999 & ...
    strcmp(aggMAE_table.trainsetSizeMax, '15*dim') & ...
    strcmp(aggMAE_table.meanFcn, 'meanConst');

% show these rows from the table
oldSettingsTable = aggMAE_table(origSettingsRows,:);
disp(oldSettingsTable);
oldSettingsHash = oldSettingsTable.hash{1};
oldID = find(strcmp(settingsHashes, oldSettingsHash));
fprintf('Old settings #:         %d\n', oldID);
fprintf('Old settings hash:      %s\n', oldSettingsHash);
fprintf('Old settings avg. rank: %.2f\n', mean(modelRanks(oldID,:)));

% show ranking of this settings
disp('Ranking of the old settings model:');
nSnapshots = length(snapshots);
dispModelRanks = zeros(nSnapshots+1, length(functions));
dispModelRanks(1,:) = functions;
for i = 1:nSnapshots
  dispModelRanks(i+1, :) = modelRanks(oldID,i:nSnapshots:end);
end
disp(dispModelRanks);
%}

