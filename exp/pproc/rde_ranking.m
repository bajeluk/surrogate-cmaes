% rde_ranking.m -- Identify the best model settings to be used in ModelPool

%% Load data
% load exp/experiments/exp_GPtest_01/modelStatistics.mat

  %% Calculate model settings ranks
  %
  dim_i = 1;

  dim = dimensions(dim_i);
  nSettings = length(folderModelOptions);
  nInstances = length(instances);

  modelRanks = zeros(nSettings, length(functions)*length(snapshots));
  nTrained   = zeros(nSettings, length(functions)*length(snapshots));
  dimId = find(dim == dimensions, 1);

  for func = functions
    for snp_i = 1:length(snapshots)
      snp = snapshots(snp_i);
      rows = (aggRDE_table.dim == dim) & (aggRDE_table.snapshot == snp);
      modelRanks(:,(func-1)*3 + snp_i) = ranking(cell2mat(aggRDE(rows, 5 + ((func-1)*3)))');
      nTrained(:,(func-1)*3 + snp_i)   = cell2mat(aggRDE(rows, 6 + ((func-1)*3)))';
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

  maxRank = 5;
  minTrainedPerc = 0.85;

  % bool vector of already covered functions/snapshots
  isCovered = false(1, length(functions)*length(snapshots));
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

  disp('The following settings were chosen (with # of hits underneath):')
  nHits = sum(boolSets(chosenSets, :), 2);
  disp([find(chosenSets)'; nHits']);
  fprintf('The total number of settings is %d.\n', sum(chosenSets));

  testedOptions = { 'covFcn', 'trainsetType', 'trainRange', 'trainsetSizeMax', 'meanFcn' };
  bestSettings = cell(sum(chosenSets), 2+length(testedOptions));
  bestSettings(:,1) = num2cell(find(chosenSets));
  bestSettings(:,2) = num2cell(nHits);
  for opi = 1:length(testedOptions)
    opt = testedOptions{opi};
    bestSettings(:, 2+opi) = cellfun(@(x) getfield(x, opt), folderModelOptions(find(chosenSets)), 'UniformOutput', false);
  end
  disp('Chosen settings:');
  disp(bestSettings);
