%% load GP model results
load('exp/experiments/exp_GPtest_01/modelStatistics.mat');

%% Core settings
dimensions = 2;
snapshots = [3, 9];
multiFieldNames = { 'covFcn', 'trainsetType', 'trainRange', 'trainsetSizeMax', 'meanFcn' };

%% Model choice settings

% maximal allowed rank for choosing the model for the function/snapshot
rdeRankingOpts.maxRank = {35, 50};
% index of the statistic by which the ranking is done, 1 = mean, 2 = 75%-quantile
rdeRankingOpts.colStat = {2};
% minimal allowed training success rate for choosing the model for the function/snapshot
rdeRankingOpts.minTrainedPerc = {0.70, 0.85};
% index of the statistic by which the error-ranking is done, 1 = mean, 2 = 75%-quantile
rdeRankingOpts.colSetCover = {1, 2};
% transformation of the number of covered functions/snapshots
% ('fsCovered' is a column vector of these number for each settings)
% and ranks of respective models for covering functions/snapshots
% ('mRanks' is a matrix of size 'nSettings x sum(~isCovered)')
rdeRankingOpts.f_weight = { @(nCover, modelErrors) (nCover.^(3))./sum(modelErrors, 2), ...
    @(nCover, modelErrors) (nCover.^(4))./sum(modelErrors, 2) };

% get the numbers of factors
optsFields = fieldnames(rdeRankingOpts);
nFactors = zeros(size(optsFields));
for f = 1:length(optsFields)
  thisField = optsFields{f};
  nFactors(f) = length(rdeRankingOpts.(thisField));
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
    opts.(thisField) = rdeRankingOpts.(thisField){ffDesign(ffi, f)};
  end
  
  % perform the "test"
  run('rde_ranking');
  bestSettingsTables{ffi} = bestSettingsTable{1};
  bestSettingsCells{ffi} = bestSettingsCell{1};
  bestSettingsStructs{ffi} = bestSettingsStruct{1};
  bestSettingsOpts{ffi} = opts;
end

%% make a table with results
rdeRankingResults_header = ['id', optsFields', ...
  { 'nSettings', 'avg_nCovered', ...
    'avg_RDE', 'avg_covrd_RDE', 'avg_rank', 'avg_covrd_rank', ...
    'avg_n_covrd_funs', 'n_ARD', 'n_Linear' }];
rdeRankingResults = table((1:n_fullfact)');
for f = 1:length(optsFields)
  rdeRankingResults(:,1+f) = cell2table( ...
       cellfun(@(x) x.(optsFields{f}), bestSettingsOpts, 'UniformOutput', false)');
end

for ffi = 1:n_fullfact
  % nSettings
  f = length(optsFields) + 2;
  rdeRankingResults{ffi, f} = size(bestSettingsTables{ffi}, 1);
  % avg_nCovered
  f = f + 1;
  rdeRankingResults{ffi, f} = nanmean(bestSettingsTables{ffi}.No_covered);
  % avg_RDE
  f = f + 1;
  rdeRankingResults{ffi, f} = nanmean(bestSettingsTables{ffi}.avg_RDE);
  % avg_covrd_RDE
  f = f + 1;
  rdeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.covrd_RDE);
  % avg_rank
  f = f + 1;
  rdeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.avg_rank);
  % avg_covrd_rank
  f = f + 1;
  rdeRankingResults{ffi, f} = mean(bestSettingsTables{ffi}.covrd_rank);
  % avg_n_covrd_funs
  f = f + 1;
  rdeRankingResults{ffi, f} = mean( cellfun( ...
      @(x) length(strsplit(x, ' ')), bestSettingsTables{ffi}.covrd_funs) );
  % n_ARD
  f = f + 1;
  rdeRankingResults{ffi, f} = sum( cellfun( @(x) ~isempty(x), ...
      regexpi(bestSettingsTables{ffi}.covFcn, 'ard') ));
  % n_Linear
  f = f + 1;
  rdeRankingResults{ffi, f} = sum( cellfun( @(x) ~isempty(x), ...
      regexpi(bestSettingsTables{ffi}.meanFcn, 'Linear') ));
end
rdeRankingResults.Properties.VariableNames = rdeRankingResults_header;

%% output results
fprintf(['\nBest sets of settings according to the average RDE\n'...
    'of covered fuctions/snapshot:\n\n']);
disp(sortrows(rdeRankingResults, 'avg_covrd_RDE'));

% print the 3 best settings options structs
N_PRINT = 3;
[~, idx] = sort(rdeRankingResults.avg_covrd_RDE);
for i = 1:N_PRINT
  fprintf('\n== Best set of settings #%d (id = %d)== \n\n', i, idx(i));
  printStructure(bestSettingsStructs{idx(i)});
end

%% find aggRDE_table rows which correspond to the historically most used settings
%{
origSettingsRows = strcmp(aggRDE_table.covFcn, '{@covMaterniso, 5}') & ...
    strcmp(aggRDE_table.trainsetType, 'clustering') & ...
    aggRDE_table.trainRange == 0.999 & ...
    strcmp(aggRDE_table.trainsetSizeMax, '15*dim') & ...
    strcmp(aggRDE_table.meanFcn, 'meanConst');

% show these rows from the table
oldSettingsTable = aggRDE_table(origSettingsRows,:);
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

