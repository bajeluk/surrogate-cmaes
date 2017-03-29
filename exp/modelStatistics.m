function [aggRDE_table, aggMSE_table, RDEs, MSEs] = modelStatistics(modelFolders, functions, dimensions, instances, snapshots, exp_id, varargin)
% modelStatistics --  creates a table with statistics,
%                     each aggRow represents model and number of snapshot (combined),
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
%   exp_id       - EXP_ID of the experiment where to store the results and of which
%                  to try to run the experiment script with modelOptions definition
%
% Name-value pair arguments:
%   'savedModelStatistics', FILENAME -- the filename from where to load RDEs, MSEs and 
%                  folderModelOptions 
%
% Output:
%   tables with average RDE and MSE values
%
% TODO:
% [ ]  make the aggregation statistics configurable (not fixed as mean,
%      std,...)
% [ ]  choose one method of saving numerical results (arrays, or cell arrays)



  % Default input parameters settings
  if (~exist('functions', 'var') || isempty(functions))
    functions = 1:24; end
  if (~exist('dimensions', 'var') || isempty(dimensions))
    dimensions = [2, 5, 10]; end
  if (~exist('instances', 'var') || isempty(instances))
    instances = [1:5, 41:50]; end
  if (~exist('snapshots', 'var') || isempty(snapshots))
    snapshots = [1 3 9]; end
  if (~exist('exp_id', 'var') || isempty(exp_id))
    exp_id = 'exp_GPtest_01'; warning('Using default exp_id = %s !', exp_id); end

  opts = settings2struct(varargin{:});
  opts.savedModelStatistics = defopts(opts, 'savedModelStatistics', []);

  assert(isnumeric(functions), '''funcToTest'' has to be integer')
  assert(isnumeric(dimensions), '''dimsToTest'' has to be integer')
  assert(isnumeric(instances), '''instToTest'' has to be integer')
  assert(isnumeric(snapshots), '''instToTest'' has to be integer')

  % take all *model_* directories in the given directory if not directories
  % *model_* given
  modelFolders = expandModelSubdirs(modelFolders);

  if (isempty(opts.savedModelStatistics))
    % load the results from individual model-testing results files

    % prepare resulting cell arrays
    RDEs = cell(length(modelFolders), length(functions), length(dimensions));
    MSEs = cell(length(modelFolders), length(functions), length(dimensions));

    isTrained = cell(length(modelFolders), length(functions), length(dimensions));
    folderModelOptions = cell(1, length(modelFolders));

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
              else
                isTrained{i_model, i_func, i_dim} = ~ isnan(MSEs{i_model, i_func, i_dim});
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

  else
    try
      fprintf(['Trying to load the RDEs, MSEs and folderModelOptions from the file\n' ...
               '  ''%s'' ...\n'], opts.savedModelStatistics);
      saved = load(opts.savedModelStatistics);
      RDEs = saved.RDEs;
      MSEs = saved.MSEs;
      folderModelOptions = saved.folderModelOptions;
      isTrained = cellfun(@(x) ~isnan(x), MSEs, 'UniformOutput', false);
    catch err
      error('ERROR: Saved statistics cannot be loaded from ''%s'': %s', ...
          opts.savedModelStatistics, err.message);
    end
  end

  rowsCount = length(modelFolders)*length(snapshots);
  % three descriptive columns for hash, dimension and snapshot
  N_DESCR_COLS = 3;
  % the number of aggreg. columns is 3-times more than # of functions,
  % because 'mean', 'std' and '# of sucess trains' are stored
  aggColumnsCount = N_DESCR_COLS + 3*length(functions);
  % the number of data columns in 'allXXX' corresponds to the number of
  % functions
  allColumnsCount = length(functions);

  % aggregated tables (aggregation via mean, standard deviation
  % and # of success trains)
  %
  % TODO: make the aggregation statistics configurable
  aggRDE    = cell(rowsCount, aggColumnsCount);
  aggMSE    = cell(rowsCount, aggColumnsCount);
  colNames  = cell(1, aggColumnsCount);
  rowNames  = cell(1, rowsCount);

  % try to run the script 'exp/experiments/[EXP_ID].m'
  % and load the multi-valued fields from the "loaded" modelOptions
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
  % use the modelOpts loaded from the first model test results above
  if (isempty(multiFieldNames))
    multiFieldNames = sort(fields(folderModelOptions{1}));
    warning('Unable to load non-factorial design of modelOptions. Using all possible fields (%d) from some modelOptions loaded from test results.', length(multiFieldNames));
  end
  % prepare cell-array for the values of the multi-valued fields
  modelOptsValues = cell(rowsCount, length(multiFieldNames));

  % non-aggregated tables (having non-aggregated all measured values)
  % the columns are:
  %     hash, dimension, snapshot#, ...[multi-valued values] ..., f1, f2, ...
  % store the values into separate numerical arrays (not cell-arrays as in
  % the case of aggregated statistics
  %
  % TODO: choose one method of saving numerical results (arrays/cell arrays)
  allRDE_cell = cell(1, N_DESCR_COLS + length(multiFieldNames));
  allRDE_mat  = NaN(1, allColumnsCount);
  allMSE_cell = cell(1, N_DESCR_COLS + length(multiFieldNames));
  allMSE_mat  = NaN(1, allColumnsCount);
  allColNames = cell(1, N_DESCR_COLS + length(multiFieldNames) + allColumnsCount);

  aggRow    = 0;
  allRow    = 0;
  N_INST    = length(instances);

  % set-up a new table aggRow for each model, snapshot and dimension
  for i_model = 1:length(modelFolders)
    [modelName, hash] = parseFolderName(modelFolders{i_model});

    for i_snapshot = 1:length(snapshots)
      snapshot = snapshots(i_snapshot);

      for i_dim = 1:length(dimensions)
        dim = dimensions(i_dim);

        aggRow = aggRow + 1;
        rowNames{aggRow} = sprintf('%s_%s_%dD_S%d', modelName, hash, dim, snapshot);

        % save values of model settings (from modelOptions) which had
        % multiple possible values
        for m_opt = 1:length(multiFieldNames)
          modelOptsValues{aggRow, m_opt} = folderModelOptions{i_model}.(multiFieldNames{m_opt});
        end

        % save the model's hash, current dimension and snapshot# into the first three columns
        colNames(1:N_DESCR_COLS) = {'hash', 'dim', 'snapshot'};
        aggRDE(aggRow, 1:N_DESCR_COLS) = {hash, dim, snapshot};
        aggMSE(aggRow, 1:N_DESCR_COLS) = {hash, dim, snapshot};

        % prepare **all values tables**
        % columns are: hash, dim, snapshot, ...[multiFieldNames]..., f1, f2, f3, ...
        allColNames(1:N_DESCR_COLS) = {'hash', 'dim', 'snapshot'};
        allColNames(N_DESCR_COLS + (1:length(multiFieldNames))) = multiFieldNames;
        allRDE_cell(allRow+(1:N_INST), 1:N_DESCR_COLS) = repmat({hash, dim, snapshot}, N_INST, 1);
        allMSE_cell(allRow+(1:N_INST), 1:N_DESCR_COLS) = repmat({hash, dim, snapshot}, N_INST, 1);
        % save modelOption-values (already here, they are the same for every
        % function (in fact even for dimension...)
        allRDE_cell(allRow+(1:N_INST), N_DESCR_COLS + (1:length(multiFieldNames))) = ...
            repmat(modelOptsValues(aggRow, :), N_INST, 1);
        allMSE_cell(allRow+(1:N_INST), N_DESCR_COLS + (1:length(multiFieldNames))) = ...
            repmat(modelOptsValues(aggRow, :), N_INST, 1);
        % prepare also matrices for numerical results (filled with NaNs)
        allRDE_mat(allRow+(1:N_INST), :) = NaN(N_INST, allColumnsCount);
        allMSE_mat(allRow+(1:N_INST), :) = NaN(N_INST, allColumnsCount);

        col = 3;
        allCol = 0;

        % save statistics from each model / function / dimension / snapshot
        for i_func = 1:length(functions)
          func = functions(i_func);

          % all values
          allCol = allCol + 1;
          allColNames{N_DESCR_COLS + length(multiFieldNames) + allCol} = sprintf('f%d', func);
          n_inst = min(N_INST, length(RDEs{i_model, i_func, i_dim}(:, i_snapshot)));
          allRDE_mat(allRow+(1:n_inst), allCol) = RDEs{i_model, i_func, i_dim}(1:n_inst, i_snapshot);
          allMSE_mat(allRow+(1:n_inst), allCol) = MSEs{i_model, i_func, i_dim}(1:n_inst, i_snapshot);

          % mean
          col = col + 1;
          colNames{col} = sprintf('f%d_m', func);
          aggRDE{aggRow, col} = nanmean(RDEs{i_model, i_func, i_dim}(:, i_snapshot));
          aggMSE{aggRow, col} = nanmean(MSEs{i_model, i_func, i_dim}(:, i_snapshot));

          % standard deviation
          col = col + 1;
          colNames{col} = sprintf('f%d_std', func);
          aggRDE{aggRow, col} = nanstd(RDEs{i_model, i_func, i_dim}(:, i_snapshot));
          aggMSE{aggRow, col} = nanstd(MSEs{i_model, i_func, i_dim}(:, i_snapshot));

          % the number of successfuly trained models
          col = col + 1;
          colNames{col} = sprintf('f%d_tr', func);
          if (~isempty(isTrained{i_model, i_func, i_dim}))
            aggRDE{aggRow, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
            aggMSE{aggRow, col} = sum(isTrained{i_model, i_func, i_dim}(:, i_snapshot));
          else
            aggRDE{aggRow, col} = sum(~isnan(MSEs{i_model, i_func, i_dim}(:, i_snapshot)));
            aggMSE{aggRow, col} = sum(~isnan(MSEs{i_model, i_func, i_dim}(:, i_snapshot)));
          end
        end
        allRow = allRow + N_INST;
      end
    end
  end

  % all values tables
  allRDE_table = cell2table(allRDE_cell);
  allRDE_table(:, end+(1:size(allRDE_mat,2))) = num2cell(allRDE_mat);
  allRDE_table.Properties.VariableNames = allColNames;

  allMSE_table = cell2table(allMSE_cell);
  allMSE_table(:, end+(1:size(allMSE_mat,2))) = num2cell(allMSE_mat);
  allMSE_table.Properties.VariableNames = allColNames;

  % aggregated tables
  aggRDE_table = cell2table([modelOptsValues, aggRDE], 'VariableNames', [multiFieldNames, colNames]);
  aggRDE_table.Properties.RowNames = rowNames;
  disp('RDE Table');
  disp(aggRDE_table);
  writetable(aggRDE_table, fullfile(exppath_short, exp_id, 'aggRDE_table.txt'), 'WriteRowNames', true);

  aggMSE_table = cell2table([modelOptsValues, aggMSE], 'VariableNames', [multiFieldNames, colNames]);
  aggMSE_table.Properties.RowNames = rowNames;
  disp('MSE Table');
  disp(aggMSE_table);
  writetable(aggMSE_table, fullfile(exppath_short, exp_id, 'aggMSE_table.txt'), 'WriteRowNames', true);

  save(fullfile(exppath_short, exp_id, 'modelStatistics.mat'), 'aggRDE_table', 'aggMSE_table', 'RDEs', 'aggRDE', 'MSEs', 'aggMSE', 'allRDE_mat', 'allRDE_table', 'allMSE_mat', 'allMSE_table', 'isTrained', 'folderModelOptions', 'modelFolders', 'functions', 'dimensions', 'snapshots', 'instances');
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

function modelFolders = expandModelSubdirs(modelFolders)
  % take all '*model_*' directories in the given directory if 'dir' is not a cell array
  % of directories  with pattern  '*model_*'
  if ((iscell(modelFolders) && length(modelFolders) == 1))
    modelFolders = modelFolders{1}; end
  if (ischar(modelFolders))
    [~, dirItself] = fileparts(modelFolders);
    if (isempty(strfind(dirItself, 'model_')))
      warning('Considering the modelFolder ''%s'' as a directory with models...', modelFolders);
      dirs = dir(fullfile(modelFolders, '*model_*'));
      dirNames = { dirs.name };
      modelFolders = cellfun(@(x) fullfile(modelFolders, x), dirNames, 'UniformOutput', false);
    else
      modelFolders = { modelFolders };
    end
  end
end
