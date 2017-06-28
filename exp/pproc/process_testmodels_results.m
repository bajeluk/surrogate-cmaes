% Parameters

expid = 'exp_DTSmodels_01';
snapshotGroups = { [5,6,7], [18,19,20] };
errorCol = 'rde2'; % 'rdeM1_M2WReplace'; % 'rdeValid';
nTrainedCol = 'nTrained2';
plotImages = true;
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

% Initialization
modelErrorRanksPerFS = cell(length(dimensions), 1);
modelErrorRanks = zeros(length(hashes), length(dimensions));
modelErrors = zeros(length(hashes), length(dimensions));
bestModelNumbers = zeros(length(hashes), length(dimensions));
bestModelRankNumbers = zeros(length(hashes), length(dimensions));
for di = 1:length(dimensions)
  modelErrorRanksPerFS{di} = zeros(length(hashes), length(functions) * nSnapshotGroups);  
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


% Summarizing results
woF5 = [1:8, 11:48];
for di = 1:length(dimensions)
  for col = woF5
    modelErrorRanksPerFS{di}(:, col) = medianRank(modelErrorsPerFS{di}(:, col));
  end
  modelErrors(:, di) = nansum(modelErrorsPerFS{di}(:, woF5), 2) ./ (length(woF5));
  modelErrorRanks(:, di) = nansum(modelErrorRanksPerFS{di}(:, woF5), 2) ./ (length(woF5));
  modelErrorsDivSuccess(:, di) = nansum( modelErrorsPerFS{di}(:, woF5) ...
      ./ trainSuccessPerFS{di}(:, woF5), 2 ) ./ (length(woF5));
  modelErrorDivSuccessRanks(:, di) = nansum( modelErrorRanksPerFS{di}(:, woF5) ...
      ./ trainSuccessPerFS{di}(:, woF5), 2 ) ./ (length(woF5));  
  % Normlize to (0, 1)
  % modelErrors(:, di) = (modelErrors(:, di) - min(modelErrors(:, di))) ./ (max(modelErrors(:, di)) - min(modelErrors(:, di)));
  [~, bestModelNumbers(:, di)] = sort(modelErrorsDivSuccess(:, di));
  [~, bestModelRankNumbers(:, di)] = sort(modelErrorDivSuccessRanks(:, di));
end  % for dimensions

if (plotImages)
  for di = 1:length(dimensions)
    dim = dimensions(di);
    figure();
    image(modelErrorsPerFS{di}(:, woF5) ./ trainSuccessPerFS{di}(:, woF5), 'CDataMapping', 'scaled');
    colorbar;
    ax = gca();
    ax.XTick = 1:2:length(woF5);
    ax.XTickLabel = ceil(woF5(1:2:end) ./ nSnapshotGroups);
    % ax.XTickLabel = cellfun(@(x) regexprep(x, '_.*', ''), modelErrorsColNames(woF5), 'UniformOutput', false);
    title([num2str(dim) 'D']);
    xlabel('functions and snapshot groups');
  end
  
  figure();
  image(modelErrorsDivSuccess, 'CDataMapping', 'scaled');
  colorbar;
  ax = gca();
  ax.XTick = 1:length(dimensions);
  ax.XTickLabel = dimensions;
  title('Average normalized RDE');
  xlabel('dimension');
end

%% Anova-n
p = {};
stats = {};
for di = 1:length(dimensions)
  if (exist('modelsSettingsFSI', 'var'))
    % Prepare Anova-n y's and categorical predictors
    categorical = { modelsSettingsFSI{di}(:, 3), cell2mat(modelsSettingsFSI{di}(:, 4)), ...
        cellfun(@(x) str2num(regexprep(x, '\*.*', '')), modelsSettingsFSI{di}(:, 5)), ...
        modelsSettingsFSI{di}(:, 6) };
    y = modelErrorsDivSuccessFSI{di};
  else
    % Prepare Anova-n y's and categorical predictors
    categorical = { modelsSettings(:, 3), cell2mat(modelsSettings(:, 4)), ...
        cellfun(@(x) str2num(regexprep(x, '\*.*', '')), modelsSettings(:, 5)), ...
        modelsSettings(:, 6) };
    y = modelErrorsDivSuccess(:, di);
  end
  % Anova-n itself:
  [p{di},tbl,stats{di},terms] = anovan(y, categorical, 'model', 1, 'varnames', multiFieldNames, 'display', 'off');
end

%% Multcompare
mstd = {};
mltc = {};
for di = 1:length(dimensions)
  fprintf('==== %d D ====\n', dimensions(di));
  for i = 1:length(multiFieldNames)
    fprintf('==   predictor: %s   ==\n', multiFieldNames{i});
    
    % Do the multi-comparison
    [c,mstd{di, i},h,nms] = multcompare(stats{di}, 'Dimension', i, 'display', 'off');
    mltc{di, i} = c;
    
    % Identify the lowest estimated mean
    [minMean, iMinMean] = min(mstd{di, i}(:, 1));
    
    % Find other values which are not statistically different
    otherLowMat = c((c(:,6) >= 0.05) & (ismember(c(:,1), iMinMean) | ismember(c(:,2), iMinMean)), 1:2);
    otherLowIdx = [];
    for j = 1:size(otherLowMat,1)
      otherLowIdx(j) = otherLowMat(j, ~ismember(otherLowMat(j,:), iMinMean));
    end
    fprintf('Lowest mean is for %s.\n', nms{iMinMean});

    if (any(otherLowIdx))
      fprintf('All lowest are:\n\n');
      disp(nms(otherLowIdx));
    else
      fprintf('No other low values.\n\n');
    end
  end
end
