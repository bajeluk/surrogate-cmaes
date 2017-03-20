function [rdeTable, mseTable, RDEs, MSEs] = modelStatistics(modelFolders, functions, dimensions, instances, snapshots)
% modelStatistics --  creates a table with statistics,
%                     each row represents model and number of snapshot (combined),
%                     each column represents function and dimension (also combined),
%                     value of cell is average of values from given instances
%
% Input:
%   modelFolders  - list of folders with models that will be displayed in table | cell-array of
%                  string
%   functions    - functions that will be displayed in table | array of integers
%   dimensions   - dimensions that will be displayed in table | array of integers
%   instances    - instances that will be displayed in table | array of integers
%   snapshots    - snapshots that will be displayed in table | array of integers
%
% Output:
%   tables with average RDE and MSE values
%

  % Default input parameters settings
  if (~exist('functions', 'var'))
    functions = 1; end
  if (~exist('dimensions', 'var'))
    dimensions = 2; end
  if (~exist('instances', 'var'))
    instances = [1:5, 41:50]; end
  if (~exist('snapshots', 'var'))
    snapshots = [1 3 9]; end

  assert(isnumeric(functions), '''funcToTest'' has to be integer')
  assert(isnumeric(dimensions), '''dimsToTest'' has to be integer')
  assert(isnumeric(instances), '''instToTest'' has to be integer')
  assert(isnumeric(snapshots), '''instToTest'' has to be integer')

  rowsCount = length(modelFolders)*length(snapshots);
  columnsCount = length(functions)*length(dimensions);

  % prepare resulting cell arrays
  RDEs = cell(length(modelFolders), length(functions), length(dimensions));
  MSEs = cell(length(modelFolders), length(functions), length(dimensions));
  isTrained = cell(length(modelFolders), length(functions), length(dimensions));

  % {
  % iterate through all specified models / dimensions / functions
  for i_model = 1:length(modelFolders)
    [modelName, hash] = parseFolderName(modelFolders{i_model});
    fprintf('Model %d [%s_%s]... ', i_model, modelName, hash);

    for i_dim = 1:length(dimensions)
      dim = dimensions(i_dim);
      fprintf('%dD: ', dim);

      for i_func = 1:length(functions)
        func = functions(i_func);
        fprintf('f%d ', func);
        is_empty = false;

        % load results for this model / fun / dim
        fileName = fullfile(modelFolders{i_model}, sprintf('%s_%s_f%d_%dD.mat', modelName, hash, func, dim));
        if (~exist(fileName, 'file'))
          warning('File %s not found', fileName);
          is_empty = true;
        else
          data = load(fileName);
          data_instances = ismember(data.instances, instances);
          if (any(data_instances))

            % save the statistics into respective resulting cell arrays
            RDEs{i_model, i_func, i_dim} = data.stats.rde(data_instances, snapshots);
            MSEs{i_model, i_func, i_dim} = data.stats.mse(data_instances, snapshots);
            isTrained{i_model, i_func, i_dim} = cellfun( ...
                @(m) m.isTrained(), data.models(data_instances, snapshots) );

          else
            warning('There''s no instance with any results for %s_%s, f%d in %dD', ...
                modelName, hash, func, dim);
            is_empty = true;
          end
        end
        if (is_empty)
          RDEs{i_model, i_func, i_dim} = [];
          MSEs{i_model, i_func, i_dim} = [];
          isTrained{i_model, i_func, i_dim} = [];
        end
      end
      fprintf('\n');
    end
  end

  % }
  % load('/tmp/stats.mat');

  valuesRDE = cell(rowsCount, columnsCount);
  valuesMSE = cell(rowsCount, columnsCount);
  colNames  = cell(1, columnsCount);
  rowNames  = cell(1, rowsCount);
  row = 0;

  % set-up a new table row for each model, snapshot and dimension
  for i_model = 1:length(modelFolders)
    [modelName, hash] = parseFolderName(modelFolders{i_model});

    for i_snapshot = 1:length(snapshots)
      snapshot = snapshots(i_snapshot);

      for i_dim = 1:length(dimensions)
        dim = dimensions(i_dim);

        row = row + 1;
        rowNames{row} = sprintf('%s_%s_%dD_S%d', modelName, hash, dim, snapshot);

        % save current dimension and snapshot into the first two columns
        colNames(1:2) = {'dim', 'snapshot'};
        valuesRDE(row, 1:2) = {dim, snapshot};
        valuesMSE(row, 1:2) = {dim, snapshot};
        col = 2;

        % save statistics from each model / function / dimension / snapshot
        for i_func = 1:length(functions)
          func = functions(i_func);

          % mean
          col = col + 1;
          colNames{col} = sprintf('f%d_m', func);
          valuesRDE{row, col} = nanmean(RDEs{i_model, i_func, i_dim}(:, i_snapshot));
          valuesMSE{row, col} = nanmean(MSEs{i_model, i_func, i_dim}(:, i_snapshot));

          % standard deviation
          col = col + 1;
          colNames{col} = sprintf('f%d_std', func);
          valuesRDE{row, col} = nanstd(RDEs{i_model, i_func, i_dim}(:, i_snapshot));
          valuesMSE{row, col} = nanstd(MSEs{i_model, i_func, i_dim}(:, i_snapshot));

          % the number of successfuly trained models
          col = col + 1;
          colNames{col} = sprintf('f%d_tr', func);
          valuesRDE{row, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
          valuesMSE{row, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
        end
      end
    end
  end

  rdeTable = cell2table(valuesRDE, 'VariableNames', colNames);
  rdeTable.Properties.RowNames = rowNames;
  disp('RDE Table');
  disp(rdeTable);

  mseTable = cell2table(valuesMSE, 'VariableNames', colNames);
  mseTable.Properties.RowNames = rowNames;
  disp('MSE Table');
  disp(mseTable);
end

function [modelName, hash, FEs] = parseFolderName(dirName)
  [~, fname] = fileparts(dirName);
  parsed = strsplit(fname, '_');
  if (length(parsed) == 3)
    modelName   = parsed{1};
    hash        = parsed{2};
    FEs         = parsed{3};
  else
    warning('Directory name cannot be parsed.');
    modelName = [];
    hash = [];
    FEs = [];
  end
end
