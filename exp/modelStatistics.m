function [rdeTable, mseTable, RDEs, MSEs] = modelStatistics(modelFolders, functions, dimensions, instances, snapshots, exp_id)
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
  if (~exist('functions', 'var') || isempty(functions))
    functions = 1; end
  if (~exist('dimensions', 'var') || isempty(dimensions))
    dimensions = 2; end
  if (~exist('instances', 'var') || isempty(instances))
    instances = [1:5, 41:50]; end
  if (~exist('snapshots', 'var') || isempty(snapshots))
    snapshots = [1 3 9]; end
  if (~exist('exp_id', 'var') || isempty(exp_id))
    exp_id = 'exp_GPtest_01'; warning('Using default exp_id = %s !', exp_id); end

  assert(isnumeric(functions), '''funcToTest'' has to be integer')
  assert(isnumeric(dimensions), '''dimsToTest'' has to be integer')
  assert(isnumeric(instances), '''instToTest'' has to be integer')
  assert(isnumeric(snapshots), '''instToTest'' has to be integer')

  % prepare resulting cell arrays
  RDEs = cell(length(modelFolders), length(functions), length(dimensions));
  MSEs = cell(length(modelFolders), length(functions), length(dimensions));
  isTrained = cell(length(modelFolders), length(functions), length(dimensions));
  folderModelOptions = cell(1, length(modelFolders));
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
          
          % load the model testing results
          data = load(fileName);
          data_instances = ismember(data.instances, instances);

          % set up correct modelOptions
          if (isempty(folderModelOptions{i_model}))
            if (iscell(data.modelOptions) && length(data.modelOptions) > 1)
              [~, dirName] = fileparts(modelFolders{i_model});
              folderModelOptions{i_model} = getThisModelOption(dirName, data.modelOptions);
            elseif (isstruct(data.modelOptions))
              folderModelOptions{i_model} = data.modelOptions;
            else
              warning('data.modelOptions are mssing or have wrong format');
            end
          end

          if (any(data_instances))

            % save the statistics into respective resulting cell arrays
            RDEs{i_model, i_func, i_dim} = data.stats.rde(data_instances, snapshots);
            MSEs{i_model, i_func, i_dim} = data.stats.mse(data_instances, snapshots);
            if (isfield(data, 'models'))
              isTrained{i_model, i_func, i_dim} = cellfun( ...
                  @(m) m.isTrained(), data.models(data_instances, snapshots) );
            end

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

  rowsCount = length(modelFolders)*length(snapshots);
  % three extra columns for hash, dimension and snapshot, and
  % 3-times more columns than functions as mean, std and # of sucess is stored
  columnsCount = 3 + (3*length(functions))*length(dimensions);

  valuesRDE = cell(rowsCount, columnsCount);
  valuesMSE = cell(rowsCount, columnsCount);
  colNames  = cell(1, columnsCount);
  rowNames  = cell(1, rowsCount);

  % try to run the script exp/experiments/exp_id.m
  % and load the multi-valued fields from modelOptions
  exppath_short = fullfile('exp', 'experiments');
  expScript = fullfile(exppath_short, [exp_id '.m']);
  multiFieldNames = {};
  if (exist(expScript, 'file'))
    run(expScript);
    if (exist('modelOptions', 'var') && isstruct(modelOptions))
      multiFieldNames = getFieldsWithMultiValues(modelOptions);
    end
  end
  % if multi-field modelOptions loading was unsuccessful,
  % use the modelOpts loaded from some model test results above
  if (isempty(multiFieldNames))
    multiFieldNames = sort(fields(folderModelOptions{1}));
    warning('Unable to load non-factorial design of modelOptions. Using all possible fields (%d) from some modelOptions loaded from test results.', length(multiFieldNames));
  end
  modelOptsValues = cell(rowsCount, length(multiFieldNames));

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

        % save values of model settings (from modelOptions) which had
        % multiple possible values
        for m_opt = 1:length(multiFieldNames)
          modelOptsValues{row, m_opt} = folderModelOptions{i_model}.(multiFieldNames{m_opt});
        end

        % save current dimension and snapshot into the first two columns
        colNames(1:3) = {'hash', 'dim', 'snapshot'};
        valuesRDE(row, 1:3) = {hash, dim, snapshot};
        valuesMSE(row, 1:3) = {hash, dim, snapshot};
        col = 3;

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
          if (~isempty(isTrained{i_model, i_func, i_dim}))
            valuesRDE{row, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
            valuesMSE{row, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
          end
        end
      end
    end
  end

  rdeTable = cell2table([modelOptsValues, valuesRDE], 'VariableNames', [multiFieldNames, colNames]);
  rdeTable.Properties.RowNames = rowNames;
  disp('RDE Table');
  disp(rdeTable);
  writetable(rdeTable, fullfile(exppath_short, exp_id, 'rdeTable.txt'), 'WriteRowNames', true);

  mseTable = cell2table([modelOptsValues, valuesMSE], 'VariableNames', [multiFieldNames, colNames]);
  mseTable.Properties.RowNames = rowNames;
  disp('MSE Table');
  disp(mseTable);
  writetable(mseTable, fullfile(exppath_short, exp_id, 'mseTable.txt'), 'WriteRowNames', true);

  save(fullfile(exppath_short, exp_id, 'modelStatistics.mat'), 'rdeTable', 'mseTable', 'RDEs', 'valuesRDE', 'MSEs', 'valuesMSE', 'isTrained', 'folderModelOptions', 'modelFolders', 'functions', 'dimensions', 'snapshots');
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

function multiFieldNames = getFieldsWithMultiValues(modelOpts)
% Get names of fields with more than one value in its cell-array fieldvalue
  assert(isstruct(modelOpts));

  fnames = fieldnames(modelOpts);
  % multiFieldNames = cell(1, length(fnames));
  multiFieldNames = {};
  for i = 1:length(fnames)
    fname = fnames{i};
    if (iscell(modelOpts.(fname)) && length(modelOpts.(fname)) > 1)
      multiFieldNames{end+1} = fname;
    end
  end
end

function [mo_struct, mo_idx] = getThisModelOption(dirName, modelOptions)
  m = regexp(dirName, '_([0-9]+)_', 'tokens');
  if (~isempty(m{1}{1}))
    dirHash = str2double(m{1}{1});

    for i = 1:length(modelOptions)
      if (str2double(modelHash(modelOptions{i})) == dirHash)
        mo_idx = i;
        mo_struct = modelOptions{i};
        return;
      end
    end
  end
  fprintf(2, 'The right index of the directory hash not found in modelOptions.\n');
  mo_idx = -1;
  mo_struct = [];
end
