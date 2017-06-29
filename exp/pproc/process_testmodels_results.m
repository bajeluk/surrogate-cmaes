% Parameters

expid = 'exp_DTSmodels_01';
snapshotGroups = { [5,6,7], [18,19,20] };
errorCol = 'rde2'; % 'rdeM1_M2WReplace'; % 'rdeValid';
nTrainedCol = 'nTrained2';
plotImages = 'mean';  % 'off', 'rank'
aggFcn = @(x) quantile(x, 0.75);

% Loading the data
if (~exist('resultTableAgg', 'var'))
  load(['exp/experiments/' expid '/modelStatistics.mat']);
end
run(['exp/experiments/' expid '.m']);

% Clear-out the non-interesting settings:
%   trainRange == 1.5 OR trainsetSizeMax == 5*dim
modelOptions.trainRange(1) = [];
modelOptions.trainsetSizeMax(1) = [];

% Processing model options
nSnapshotGroups = length(snapshotGroups);
multiFieldNames = getFieldsWithMultiValues(modelOptions);
modelOptions_fullfact = combineFieldValues(modelOptions);
hashes = cellfun(@(x) modelHash(x), modelOptions_fullfact, 'UniformOutput', false);

% Model settings differences
modelsSettings = cell(length(hashes), length(multiFieldNames));
for mi = 1:length(hashes)
  modelsSettings(mi, 1:2) = {mi, hashes{mi}};
  for mf = 1:length(multiFieldNames)
    modelsSettings{mi, 2+mf} = modelOptions_fullfact{mi}.(multiFieldNames{mf});
  end
end

% Column names
modelErrorsColNames = cell(length(functions)*nSnapshotGroups, 1);
for fi = 1:length(functions)
  fun = functions(fi);
  for si = 1:nSnapshotGroups
    iCol = (fi-1)*nSnapshotGroups + si;
    modelErrorsColNames{iCol} = ['f' num2str(fun) '_S' num2str(si)];
  end
end

% Loading errors from the table
if (~exist('modelErrorsPerFS', 'var') || ~exist('modelsSettingsFSI', 'var'))
  modelErrorsPerFS = cell(length(dimensions), 1);

  for di = 1:length(dimensions)
    dim = dimensions(di);
    fprintf('%dD: ', dim);
    modelErrorsPerFS{di} = zeros(length(hashes), length(functions) * nSnapshotGroups);
    modelErrorsDivSuccessFSI{di} = [];
    modelsSettingsFSI{di} = cell(0, size(modelsSettings, 2));
    funDimSngFSI{di}   = zeros(0, 3);
    rowFSI = 0;

    for mi = 1:length(hashes)
      hash = hashes{mi};
      fprintf('model %s ', hash);

      for fi = 1:length(functions)
        fun = functions(fi);
        % fprintf('f%d ', fun);

        for si = 1:nSnapshotGroups

          iCol = (fi-1)*nSnapshotGroups + si;
          % modelErrorsPerFS{di}(mi, iCol) = resultTableAgg{ ...
          %     strcmpi(resultTableAgg.hash, hash) ...
          %     &  resultTableAgg.dim == dim ...
          %     &  resultTableAgg.fun == fun ...
          %     &  resultTableAgg.snpGroup == si, errorCol};
          thisErrors = resultTableAll{ ...
              strcmpi(resultTableAll.hash, hash) ...
              &  ismember(resultTableAll.snapshot, snapshotGroups{si}) ...
              &  resultTableAll.dim == dim ...
              &  resultTableAll.fun == fun, errorCol};
          nErrors = sum(~isnan(thisErrors));

          modelErrorsPerFS{di}(mi, iCol) = aggFcn(thisErrors);
          modelErrorsDivSuccessFSI{di}(rowFSI + (1:nErrors), 1) = ...
              thisErrors(~isnan(thisErrors)) ...
              * ((length(snapshotGroups{si})*length(instances)) / nErrors);
          [~, thisSettingsIdx] = ismember(hash, hashes);
          modelsSettingsFSI{di}(rowFSI + (1:nErrors), :) = ...
              repmat(modelsSettings(thisSettingsIdx, :), nErrors, 1);
          funDimSngFSI{di}(rowFSI + (1:nErrors), :) = ...
              repmat([fun, dim, si], nErrors, 1);
          rowFSI = rowFSI + nErrors;
        end  % for snapshotGroups

      end  % for functions

      fprintf('\n');
    end  % for models
  end
end

% Loading # of successful trains from the table
if (~exist('trainSuccessPerFS', 'var'))
  trainSuccessPerFS = cell(length(dimensions), 1);

  for di = 1:length(dimensions)
    dim = dimensions(di);
    fprintf('%dD: ', dim);
    trainSuccessPerFS{di} = zeros(length(hashes), length(functions) * nSnapshotGroups);

    for mi = 1:length(hashes)
      hash = hashes{mi};
      fprintf('model %s ', hash);

      for fi = 1:length(functions)
        fun = functions(fi);
        % fprintf('f%d ', fun);

        for si = 1:nSnapshotGroups

          iCol = (fi-1)*nSnapshotGroups + si;
          trainSuccessPerFS{di}(mi, iCol) = resultTableAgg{ ...
              strcmpi(resultTableAgg.hash, hash) ...
              &  resultTableAgg.dim == dim ...
              &  resultTableAgg.fun == fun ...
              &  resultTableAgg.snpGroup == si, nTrainedCol} ...
              / (length(snapshotGroups{si})*length(instances));
        end  % for snapshotGroups
      end  % for functions
      fprintf('\n');
    end  % for models
  end
end


%% Summarizing results
% Initialization
modelErrorRanks = zeros(length(hashes), length(dimensions)*nSnapshotGroups);
modelErrors     = zeros(length(hashes), length(dimensions*nSnapshotGroups));
bestModelNumbers = zeros(length(hashes), length(dimensions)*nSnapshotGroups);
bestModelRankNumbers = zeros(length(hashes), length(dimensions)*nSnapshotGroups);
modelErrorsDivSuccess = zeros(length(hashes), length(dimensions)*nSnapshotGroups);
modelErrorDivSuccessRanks = zeros(length(hashes), length(dimensions)*nSnapshotGroups);
modelErrorRanksPerFS = cell(length(dimensions), 1);
for di = 1:length(dimensions)
  modelErrorRanksPerFS{di} = zeros(length(hashes), length(functions) * nSnapshotGroups);  
end

for di = 1:length(dimensions)
  for si = 1:nSnapshotGroups
    idx = ((di - 1) * nSnapshotGroups) + si;

    woF5 = ((setdiff(functions, 5) - 1) * nSnapshotGroups) + si;
    for col = woF5
      modelErrorRanksPerFS{di}(:, col) = medianRank(modelErrorsPerFS{di}(:, col));
    end
    modelErrors(:, idx) = nansum(modelErrorsPerFS{di}(:, woF5), 2) ./ (length(woF5));
    modelErrorRanks(:, idx) = nansum(modelErrorRanksPerFS{di}(:, woF5), 2) ./ (length(woF5));
    modelErrorsDivSuccess(:, idx) = nansum( modelErrorsPerFS{di}(:, woF5) ...
        ./ trainSuccessPerFS{di}(:, woF5), 2 ) ./ (length(woF5));
    modelErrorDivSuccessRanks(:, idx) = nansum( modelErrorRanksPerFS{di}(:, woF5) ...
        ./ trainSuccessPerFS{di}(:, woF5), 2 ) ./ (length(woF5));  
    % Normlize to (0, 1)
    % modelErrors(:, idx) = (modelErrors(:, idx) - min(modelErrors(:, idx))) ./ (max(modelErrors(:, idx)) - min(modelErrors(:, idx)));
    [~, bestModelNumbers(:, idx)] = sort(modelErrorsDivSuccess(:, idx));
    [~, bestModelRankNumbers(:, idx)] = sort(modelErrorDivSuccessRanks(:, idx));
  end
end  % for dimensions


%% Plot images
woF5 = repelem((setdiff(functions, 5) - 1) * nSnapshotGroups, 1, nSnapshotGroups) ...
       +  repmat(1:nSnapshotGroups, 1, length(setdiff(functions, 5)));
if (~strcmpi(plotImages, 'off'))
  for di = 1:length(dimensions)
    dim = dimensions(di);
    figure();
    switch lower(plotImages)
      case 'mean'
        image(modelErrorsPerFS{di}(:, woF5) ./ trainSuccessPerFS{di}(:, woF5), 'CDataMapping', 'scaled');
      case 'rank'
        image(modelErrorRanksPerFS{di}(:, woF5) ./ trainSuccessPerFS{di}(:, woF5), 'CDataMapping', 'scaled');
      otherwise
        warning('plot style not known: %s', plotImages)
    end
    colorbar;
    ax = gca();
    ax.XTick = 1:nSnapshotGroups:size(modelErrorsPerFS{di}, 2);
    ax.XTickLabel = setdiff(functions, 5); % ceil(woF5(1:2:end) ./ nSnapshotGroups);
    % ax.XTickLabel = cellfun(@(x) regexprep(x, '_.*', ''), modelErrorsColNames(woF5), 'UniformOutput', false);
    title([num2str(dim) 'D']);
    xlabel('functions and snapshot groups');
  end
  
  figure();
  switch lower(plotImages)
    case 'mean'
      image(modelErrorsDivSuccess, 'CDataMapping', 'scaled');
    case 'rank'
      image(modelErrorDivSuccessRanks, 'CDataMapping', 'scaled');
  end
  colorbar;
  ax = gca();
  ax.XTick = 1:2:size(modelErrorsDivSuccess, 2);
  ax.XTickLabel = dimensions;
  title('RDE averaged across functions');
  xlabel('dimension, snapshotGroup');
end

%% Anova-n
p = {};
stats = {};
snapshotGroupCol = 3;
for di = 1:length(dimensions)
  for si = 1:nSnapshotGroups
    idx = ((di - 1) * nSnapshotGroups) + si;
    if (exist('modelsSettingsFSI', 'var'))
      % Prepare Anova-n y's and categorical predictors
      thisSNG = funDimSngFSI{di}(:,snapshotGroupCol) == si;
      categorical = { modelsSettingsFSI{di}(thisSNG, 3), cell2mat(modelsSettingsFSI{di}(thisSNG, 4)), ...
          cellfun(@(x) str2num(regexprep(x, '\*.*', '')), modelsSettingsFSI{di}(thisSNG, 5)), ...
          modelsSettingsFSI{di}(thisSNG, 6) };
      y = modelErrorsDivSuccessFSI{di}(thisSNG, 1);
    else
      % Prepare Anova-n y's and categorical predictors
      categorical = { modelsSettings(:, 3), cell2mat(modelsSettings(:, 4)), ...
          cellfun(@(x) str2num(regexprep(x, '\*.*', '')), modelsSettings(:, 5)), ...
          modelsSettings(:, 6) };
      y = modelErrorsDivSuccess(:, idx);
    end
    % Anova-n itself:
    [p{idx},tbl,stats{idx},terms] = anovan(y, categorical, 'model', 1, 'varnames', multiFieldNames, 'display', 'off');
  end
end

%% Multcompare
cellBestValues = cell(0, 2+length(multiFieldNames));
mstd = {};
mltc = {};
for di = 1:length(dimensions)
  for si = 1:nSnapshotGroups
    idx = ((di - 1) * nSnapshotGroups) + si;

    fprintf('==== %d D, snapshotG: %d ====\n', dimensions(di), si);
    rowStart = size(cellBestValues, 1);

    for i = 1:length(multiFieldNames)
      fprintf('==   predictor: %s   ==\n', multiFieldNames{i});

      % Do the multi-comparison
      [c,mstd{idx, i},h,nms] = multcompare(stats{idx}, 'Dimension', i, 'display', 'off');
      nValues = size(nms, 1);
      mltc{idx, i} = c;

      % Identify the lowest estimated mean (and sort them, too)
      [~, sortedMeansId] = sort(mstd{idx, i}(:, 1));
      iMinMean = sortedMeansId(1);
      fprintf('Lowest mean f-value for this predictor is for %s.\n', nms{iMinMean});

      % Find other values which are not statistically different from the
      % lowest
      otherLowRows = (c(:,6) >= 0.05) & (ismember(c(:,1), iMinMean) | ismember(c(:,2), iMinMean));
      otherLowMat = c(otherLowRows, 1:2);
      otherLowBool = false(nValues, 1);
      for j = 1:size(otherLowMat,1)
        otherLowBool(otherLowMat(j, ~ismember(otherLowMat(j,:), iMinMean))) = true;
      end
      otherLowIdx = find(otherLowBool);
      isOtherInSorted = ismember(sortedMeansId, otherLowIdx);
      otherLowSorted = sortedMeansId(isOtherInSorted);

      bestValues = cellfun(@(x) regexprep(x, '^.*=', ''), ...
          nms([iMinMean; otherLowSorted]), 'UniformOutput', false);
      cellBestValues(rowStart + (1:length(bestValues)), 2+i) = bestValues;

      if (any(otherLowBool))
        fprintf('Other statistically also low are:\n\n');
        disp(nms(otherLowSorted));
      else
        fprintf('No other statistically lowest values.\n\n');
      end
    end
    cellBestValues(end+1, :) = [{[], []}, num2cell(p{idx}')];
    cellBestValues((rowStart+1):end, 1:2) = repmat({dimensions(di), si}, ...
        size(cellBestValues,1)-rowStart, 1);
  end
end
tBestValues = cell2table(cellBestValues, ...
    'VariableNames', [{'dim', 'snG' }, multiFieldNames]);
