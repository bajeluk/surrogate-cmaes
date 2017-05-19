% rde_ranking_MP.m -- Identify the best model settings from ModelPool models

%% Load data
% load exp/experiments/exp_MPtest_01_rde/modelStatistics.mat
aggRDE_nHeaderCols = 3;
aggRDE_nColsPerFunction = 3;
nInstances = length(instances);
nSettings = length(folderModelOptions);
modelMeanRDE = zeros(nSettings, length(dimensions));
modelQuantileRDE = zeros(nSettings, length(dimensions));

for dim_i = 1:length(dimensions)

  dim = dimensions(dim_i);
  settingsHashes = cellfun(@modelHash, folderModelOptions, 'UniformOutput', false)';

  for m = 1:nSettings
    hash = settingsHashes{m};
    rows = (aggRDE_table.dim == dim) & cellfun(@(x) strcmpi(x, hash), aggRDE_table.hash);
    rdeMatrix = cell2mat(aggRDE(rows, (aggRDE_nHeaderCols+1):aggRDE_nColsPerFunction:end));
    modelMeanRDE(m, dim_i) = mean(mean(rdeMatrix));
  end
  for m = 1:nSettings
    hash = settingsHashes{m};
    rows = (aggRDE_table.dim == dim) & cellfun(@(x) strcmpi(x, hash), aggRDE_table.hash);
    rdeMatrix = cell2mat(aggRDE(rows, (aggRDE_nHeaderCols+2):aggRDE_nColsPerFunction:end));
    modelQuantileRDE(m, dim_i) = mean(mean(rdeMatrix));
  end

  settingsIndex = (1:nSettings)';
  modelOptions = cell(nSettings,5);
  for m = 1:nSettings
    modelOptions{m,1} = settingsHashes{m};
    modelOptions{m,2} = folderModelOptions{1,m}.bestModelSelection;
    modelOptions{m,3} = folderModelOptions{1,m}.historyLength;
    modelOptions{m,4} = folderModelOptions{1,m}.minTrainedModelsPercentilForModelChoice;
    modelOptions{m,5} = folderModelOptions{1,m}.maxGenerationShiftForModelChoice;
  end
  modelOptionsTable = cell2table(modelOptions);
  modelOptionsTable.Properties.VariableNames = {'hash','bestModelSelection', 'historyLength', 'minTrainedModelsPercentilForModelChoice', 'maxGenerationShiftForModelChoice'};

  resultsTable = [table(settingsIndex,modelMeanRDE, modelQuantileRDE), modelOptionsTable];

  nBestModels = 3;

  resultsTable = sortrows(resultsTable, 'modelMeanRDE');
  fprintf('%d best models by RDE mean:\n',nBestModels);
  disp(resultsTable(1:nBestModels,:));

  resultsTable = sortrows(resultsTable, 'modelQuantileRDE');
  fprintf('%d best models by RDE quantile:\n',nBestModels);
  disp(resultsTable(1:nBestModels,:));

end
