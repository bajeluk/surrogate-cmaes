% rde_ranking_MP.m -- Identify the best model settings from ModelPool models

%% Load data
% load exp/experiments/exp_MPtest_01_rde/modelStatistics.mat
aggMSE_nHeaderCols = 3;
aggMSE_nColsPerFunction = 3;
nInstances = length(instances);
nSettings = length(folderModelOptions);
modelMeanMSE = zeros(nSettings, length(dimensions));
modelQuantileMSE = zeros(nSettings, length(dimensions));

for dim_i = 1:length(dimensions)

  dim = dimensions(dim_i);
  settingsHashes = cellfun(@modelHash, folderModelOptions, 'UniformOutput', false)';

  for m = 1:nSettings
    hash = settingsHashes{m};
    rows = (aggMSE_table.dim == dim) & cellfun(@(x) strcmpi(x, hash), aggMSE_table.hash);
    mseMatrix = cell2mat(aggMSE(rows, (aggMSE_nHeaderCols+1):aggMSE_nColsPerFunction:end));
    modelMeanMSE(m, dim_i) = mean(mean(mseMatrix));
  end
  for m = 1:nSettings
    hash = settingsHashes{m};
    rows = (aggMSE_table.dim == dim) & cellfun(@(x) strcmpi(x, hash), aggMSE_table.hash);
    mseMatrix = cell2mat(aggMSE(rows, (aggMSE_nHeaderCols+2):aggMSE_nColsPerFunction:end));
    modelQuantileMSE(m, dim_i) = mean(mean(mseMatrix));
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

  resultsTable = [table(settingsIndex,modelMeanMSE, modelQuantileMSE), modelOptionsTable];

  nBestModels = 3;

  resultsTable = sortrows(resultsTable, 'modelMeanMSE');
  fprintf('%d best models by MSE mean:\n',nBestModels);
  disp(resultsTable(1:nBestModels,:));

  resultsTable = sortrows(resultsTable, 'modelQuantileMSE');
  fprintf('%d best models by MSE quantile:\n',nBestModels);
  disp(resultsTable(1:nBestModels,:));

end
