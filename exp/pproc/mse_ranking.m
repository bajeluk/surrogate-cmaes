% mse_ranking.m -- Identify the best model settings to be used in ModelPool

%% Load data
% load exp/experiments/exp_GPtest_01/modelStatistics.mat

defaultParameterSets = struct( ...
  'trainAlgorithm', { {'fmincon'} }, ...
  'hyp',            { {struct('lik', log(0.01), 'cov', log([0.5; 2]))} });
printBestSettingsDefinition = false;

maxRank = 20;
minTrainedPerc = 0.9;

bestSettings = cell(length(dimensions),1);
aggMSE_nHeaderCols = 3;
aggMSE_nColsPerFunction = 3;
nSettings = length(folderModelOptions);
modelMeanMSE = zeros(nSettings, length(dimensions));
modelMeanMSECovered = zeros(nSettings, length(dimensions));
modelMeanRankCovered = zeros(nSettings, length(dimensions));
modelFunctionsCovered = cell(nSettings, length(dimensions));

for dim_i = 1:length(dimensions)

  %% Calculate model settings ranks
  %
  dim = dimensions(dim_i);
  nInstances = length(instances);

  nSnapshots = length(snapshots);
  modelRanks = zeros(nSettings, length(functions)*nSnapshots);
  nTrained   = zeros(nSettings, length(functions)*nSnapshots);
  dimId = find(dim == dimensions, 1);
  functions(functions>100) = functions(functions>100) - 100; % for noisy functions

  for func = functions
    for snp_i = 1:nSnapshots
      snp = snapshots(snp_i);
      % take rows which belong to the considered dimension/snapshot
      % combination
      rows = (aggMSE_table.dim == dim) & (aggMSE_table.snapshot == snp);
      % take ranks according to the 2nd statistic
      colStat = 2;
      modelRanks(:,(func-1)*nSnapshots + snp_i) = ranking(cell2mat(aggMSE(rows, aggMSE_nHeaderCols+colStat + ((func-1)*aggMSE_nColsPerFunction)))');
      % # of train success should be in the 3rd column
      colNTrained = 3;
      nTrained(:,(func-1)*nSnapshots + snp_i)   = cell2mat(aggMSE(rows, aggMSE_nHeaderCols+colNTrained + ((func-1)*aggMSE_nColsPerFunction)))';
    end
  end

  %% Calculate the number of ranks 1, 2, 3,... across all functions and
  %  snapshots for each model setting

  nBestRanks = 5;

  modelBestSettings = zeros(nSettings, 1+nBestRanks);
  modelBestSettings(:,1) = 1:nSettings;
  for m = 1:nSettings
    for rnk = 1:nBestRanks
      modelBestSettings(m, 1+rnk) = sum(modelRanks(m, :) == rnk);
    end
  end

  % %% First try of model choose
  %
  % sortedSettings = sortrows(modelBestSettings, 2:6);
  %
  % % calculate points: 3p for #1, 2p for #2, 1p for #3:
  % points = bsxfun(@times, sortedSettings, [1 5 4 3 2 1]);
  % % sum the points for each model setting:
  % sumPoints = [points(:,1) sum(points(:,2:end), 2)];
  % [~, sumPoints_sorti] = sort(sumPoints(:,2), 'descend');
  %
  % nBest = 36;
  % bestModelRanks = modelRanks(sumPoints(sumPoints_sorti(1:nBest), 1), :);
  % bestModelRanks(bestModelRanks > 20) = 21;
  % bestModelRanks = 21 - bestModelRanks;
  %
  % % disp(bestModelRanks);
  % % disp(sum(bestModelRanks, 1));

  %% Identify settings wich have at least one rank #1--#5
  isRank1to5 = find(any(modelBestSettings(:,2:end) > 0, 2));

  %% Set Cover problem
  %  Identify the least number of model settings for modelPool such that
  %  these models performed (very) well on at least one function/snapshot
  %
  %  All function and snapshots has to be covered with at least one setting
  %  which performed (very) well, particularly with maximal rank 'maxRank'
  %  and the success rate of training has to be at least 'minTrainedPerc'

  % bool vector of already covered functions/snapshots
  isCovered = false(1, length(functions)*nSnapshots);
  % bool array indicating whether model settings performed (very) well on
  % respective function/snapshot
  boolSets = (modelRanks <= maxRank) & (nTrained/nInstances >= minTrainedPerc);
  % these function/snaphots can be covered by at least one settings
  possibleCover = any(boolSets, 1);
  % so-far chosen settings
  chosenSets = false(nSettings, 1);
  % sort setting-numbers according to the number of ranks #1,#2,...,#5
  sortedSettings = sortrows(modelBestSettings, 2:6);
  sortedSettings = sortedSettings(end:-1:1, :);
  sortedSettingsIdx = sortedSettings(:,1);

  while (any(~isCovered(possibleCover)))
    % identify how many uncovered sets the setting covers
    howMuchWillCover = sum(boolSets(:, ~isCovered), 2);
    % sets the number to zero for already chosen sets
    howMuchWillCover(chosenSets) = 0;
    % take all the settings with the maximal covering property
    nMaxCovers = max(howMuchWillCover);
    maxCovered = find(howMuchWillCover == nMaxCovers);
    % choose one of the max-covering sets:
    %   take such settings which has maximal covering number
    %   and is the first in sortedSettings (see above)
    [~, maxCoveredInSorted] = ismember(maxCovered, sortedSettingsIdx);
    [~, bestMaxCoveredIdx]  = min(maxCoveredInSorted);
    maxCoveringSet = maxCovered(bestMaxCoveredIdx);
    chosenSets(maxCoveringSet) = true;
    % fprintf('Picking up set #%d\n', newlyChosen);
    % set the chosen set as chosen ;)
    isCovered(boolSets(maxCoveringSet, :)) = true;
  end

  fprintf('The following settings were chosen for %dD (with # of hits underneath):', dim)
  nHits = sum(boolSets(chosenSets, :), 2);
  [find(chosenSets)'; nHits']
  fprintf('The total number of settings in %dD is %d.\n', dim, sum(chosenSets));

  % Calculate statistics of each model
  % take mean MSE statistic into modelMeanMSE (mean is calculated from instances)
  colMean = 1;
  settingsHashes = cellfun(@modelHash, folderModelOptions, 'UniformOutput', false)';

  for m = 1:nSettings
    hash = settingsHashes{m};

    % take all rows corresponding to the actual settings' hash
    rows = (aggMSE_table.dim == dim) & cellfun(@(x) strcmpi(x, hash), aggMSE_table.hash);
    mseMatrix = cell2mat(aggMSE(rows, (aggMSE_nHeaderCols+colMean):aggMSE_nColsPerFunction:end));
    modelMeanMSE(m, dim_i) = mean(mean(mseMatrix));

    % calculate also statistics based only on the functions/snapshots
    % which are covered by respective settings
    mseMatrixCovered = mseMatrix;
    modelFunctionsCovered{m, dim_i} = false(1,length(functions));
    for sni = 1:nSnapshots
      mseMatrixCovered(sni, ~boolSets(m, sni:nSnapshots:end)) = NaN;
      modelFunctionsCovered{m, dim_i} = modelFunctionsCovered{m, dim_i} ...
          | boolSets(m, sni:nSnapshots:end);
    end
    modelMeanMSECovered(m, dim_i) = mean(mseMatrixCovered(~isnan(mseMatrixCovered)));
    modelMeanRankCovered(m, dim_i) = mean(modelRanks(m, boolSets(m, :)));
  end

  headerCols    = { 'settingNo', 'No_covered', 'avg_MSE', 'covrd_MSE', 'avg_rank', 'covrd_rank', 'covrd_funs' };
  % multiFieldNames = { 'covFcn', 'trainsetType', 'trainRange', 'trainsetSizeMax', 'meanFcn' };

  nHeaderCols = length(headerCols);
  bestSettings{dim_i} = cell(1+sum(chosenSets), nHeaderCols + length(multiFieldNames));
  bestSettings{dim_i}(1, :) = [headerCols multiFieldNames];
  bestSettings{dim_i}(2:end,1) = num2cell(find(chosenSets));        % indices of chosen models
  bestSettings{dim_i}(2:end,2) = num2cell(nHits);   % the numbers of hits of each model
  bestSettings{dim_i}(2:end,3) = num2cell(modelMeanMSE(chosenSets, dim_i));       % mean of MSE accros f/snp
  bestSettings{dim_i}(2:end,4) = num2cell(modelMeanMSECovered(chosenSets, dim_i));       % mean of MSE accros f/snp which are covered by this settings
  bestSettings{dim_i}(2:end,5) = num2cell(mean(modelRanks(chosenSets, :), 2));      % mean ranks
  bestSettings{dim_i}(2:end,6) = num2cell(modelMeanRankCovered(chosenSets, dim_i));   % mean ranks of the covered f/snp
  bestSettings{dim_i}(2:end,7) = cellfun(@(x) num2str(find(x)), ...
      modelFunctionsCovered(chosenSets, dim_i), 'UniformOutput', false);  % list of covered functions

  for opi = 1:length(multiFieldNames)
    opt = multiFieldNames{opi};
    bestSettings{dim_i}(2:end, nHeaderCols + opi) = cellfun(@(x) getfield(x, opt), folderModelOptions(chosenSets), 'UniformOutput', false);
  end
  disp('Chosen settings:');
  disp(bestSettings{dim_i});

  % extract the best settings itself in a separate new cell
  % and add also the standard options from 'defaultParameterSets' defined
  % at the beginning of this file
  bestSettingsCell{dim_i} = bestSettings{dim_i}(2:end, (nHeaderCols+1):end);
  for fld = fieldnames(defaultParameterSets)'
    value = defaultParameterSets.(fld{1});
    bestSettingsCell{dim_i}(:, end+1) = value( ones(1, sum(chosenSets)), 1 );
  end

  % convert the created cell into a structure
  % and print it on the screen
  bestSettingsStruct{dim_i} = cell2struct(bestSettingsCell{dim_i}, ...
      [multiFieldNames, fieldnames(defaultParameterSets)'], 2);
  if (exist('printBestSettingsDefinition', 'var') && islogical(printBestSettingsDefinition) ...
      && printBestSettingsDefinition)
    printStructure(bestSettingsStruct{dim_i}, 'StructName', 'parameterSets');
  else
    disp('Printing the definitions of the best settings structs is switched off.')
  end

  bestSettingsTable{dim_i} = cell2table(bestSettings{dim_i}(2:end, :), ...
      'VariableNames', bestSettings{dim_i}(1,:));
  fprintf('\n-----------------------------------\n');
end
