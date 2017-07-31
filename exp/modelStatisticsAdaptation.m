function [resultTableAll, resultTableAgg] = modelStatisticsAdaptation(modelFolders, functions, dimensions, instances, snapshots, exp_id, varargin)
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
%   'firstStatistic'  -- the first statistic to calculate from set of
%                        results on (usually 15) instances (default: @mean)
%   'firstStatisticName' -- the title to be printed in table headers
%   'secondStatistic'  -- the first statistic to calculate from set of
%                        results on (usually 15) instances (default: @std)
%   'secondStatisticName' -- the title to be printed in table headers
%
% Output:
%   tables with average RDE and MSE values, and w/o aggregation
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

  if (~isempty(varargin) && ~isempty(varargin{1}) && isstruct(varargin{1}))
    opts = varargin{1};
  else
    opts = settings2struct(varargin{:});
  end
  opts.savedModelStatistics = defopts(opts, 'savedModelStatistics', []);
  opts.statistics = defopts(opts, 'statistics', { 'mse', 'rde', 'rdeValid'});
  opts.aggFunction = defopts(opts, 'aggFunction', @nanmean);
  % group of snapshots which should be aggregated together
  opts.aggSnapshots = defopts(opts, 'aggSnapshots', {snapshots});

  % opts.firstStatistic = defopts(opts, 'firstStatistic', @mean);
  % opts.firstStatisticName = defopts(opts, 'firstStatisticName', func2str(opts.firstStatistic));
  % opts.secondStatistic = defopts(opts, 'secondStatistic', @std);
  % opts.secondStatisticName = defopts(opts, 'secondStatisticName', func2str(opts.secondStatistic));

  % try to run the script 'exp/experiments/[EXP_ID].m'
  % and load the multi-valued fields from the "loaded" modelOptions
  exppath_short = fullfile('exp', 'experiments');
  expScript = fullfile(exppath_short, [exp_id '.m']);
  if (exist(expScript, 'file'))
    run(expScript);
  end

  assert(isnumeric(functions), '''funcToTest'' has to be integer')
  assert(isnumeric(dimensions), '''dimsToTest'' has to be integer')
  assert(isnumeric(instances), '''instToTest'' has to be integer')
  assert(isnumeric(snapshots), '''snpsToTest'' has to be integer')

  % take all *model_* directories in the given directory if not directories
  % *model_* given
  modelFolders = expandModelSubdirs(modelFolders);
  if (isempty(opts.savedModelStatistics))
    % load the results from individual model-testing results files

    % prepare resulting struct with cell arrays
    results = struct();
    for s = opts.statistics
      if (iscell(s) && ~isempty(s))
        fieldName = s{1};
      else
        continue;
      end
      results.(fieldName) = {};
    end

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

          % load results for this model / func / dim
          fileName = fullfile(modelFolders{i_model}, sprintf('%s_%s_f%d_%dD.mat', modelName, hash, func, dim));
          if (~exist(fileName, 'file'))
            warning('File %s not found', fileName);
            is_empty = true;
            continue;
          else

            % load the model testing results
            data = load(fileName);
            data_instances = ismember(data.instances, instances);
            max_instance = min(size(data.stats.mse, 1), find(~all(isnan(data.stats.mse), 2), 1, 'last'));
            data_instances((max_instance+1):end) = false;

            % set up correct modelOptions
            if (isempty(folderModelOptions{i_model}))
              if (iscell(data.modelOptions) && length(data.modelOptions) > 1)
                [~, dirName] = fileparts(modelFolders{i_model});
                folderModelOptions{i_model} = getThisModelOption(dirName, data.modelOptions);
              elseif (isstruct(data.modelOptions))
                folderModelOptions{i_model} = data.modelOptions;
              elseif (isstruct(data.modelOptions{1}))
                folderModelOptions{i_model} = data.modelOptions{1};
              else
                warning('data.modelOptions are mssing or have wrong format');
              end
            end

            if (any(data_instances))

              % save the statistics into resulting cell arrays
              for s = opts.statistics
                if (iscell(s) && ~isempty(s))
                  fieldName = s{1};
                else
                  continue;
                end
                results.(fieldName){i_model, i_func, i_dim} = data.stats.(fieldName)(data_instances, snapshots);
              end

              if (isfield(data, 'models'))
                results.isTrained{i_model, i_func, i_dim} = cellfun( ...
                    @(m) ~isempty(m) && m.isTrained(), data.models(data_instances, snapshots) );
              else
                results.isTrained{i_model, i_func, i_dim} = ~isnan(results.mse{i_model, i_func, i_dim});
              end
              if (isfield(data, 'models2'))
                results.isTrained2{i_model, i_func, i_dim} = cellfun( ...
                    @(m) ~isempty(m) && m.isTrained(), data.models2(data_instances, snapshots) );
              else
                results.isTrained2{i_model, i_func, i_dim} = isfield(results, 'rde2') && ~isnan(results.rde2{i_model, i_func, i_dim});
              end

            else
              warning('There''s no instance with any results for %s_%s, f%d in %dD', ...
                  modelName, hash, func, dim);
              is_empty = true;
            end
          end  % else/if (~exist(fileName, 'file'))
          if (is_empty)
            for s = opts.statistics
              results.(s{1}){i_model, i_func, i_dim} = [];
            end
            results.isTrained{i_model, i_func, i_dim} = [];
            results.isTrained2{i_model, i_func, i_dim} = [];
          end
        end
        fprintf('\n');
      end
    end

    save(fullfile(exppath_short, exp_id, ['modelResults_' datestr(now(), 'yyyy-mm-dd') '.mat']), 'results', 'folderModelOptions', 'modelFolders', 'functions', 'dimensions', 'snapshots', 'instances');

  else  % if (isempty(opts.savedModelStatistics))
    try
      fprintf(['Trying to load the model training/testing statistics from the file\n' ...
               '  ''%s'' ...\n'], opts.savedModelStatistics);
      saved = load(opts.savedModelStatistics);
      results = saved.results;
      folderModelOptions = saved.folderModelOptions;
    catch err
      error('ERROR: Saved statistics cannot be loaded from ''%s'': %s', ...
          opts.savedModelStatistics, err.message);
    end
  end

  descrColNamesAll = {'hash', 'dim', 'fun', 'inst', 'snapshot'};
  descrColNamesAgg = {'hash', 'dim', 'fun', 'snpGroup'};

  resultCellAll = cell(1, length(descrColNamesAll));
  resultMatAll  = NaN(1, length(opts.statistics) + 2);
  resultCellAgg = cell(1, length(descrColNamesAgg));
  resultMatAgg  = NaN(1, length(opts.statistics) + 2);
  resultColNamesAll = [descrColNamesAll, opts.statistics, {'isTrained', 'isTrained2'}];
  resultColNamesAgg = [descrColNamesAgg, opts.statistics, {'nTrained', 'nTrained2'}];
  resultRowNamesAll = {};
  resultRowNamesAgg = {};

  allRow = 0;
  lastFilledRow = 0;
  aggRow = 0;

  for i_model = 1:length(modelFolders)
    [modelName, hash] = parseFolderName(modelFolders{i_model});

    for i_dim = 1:length(dimensions)
      dim = dimensions(i_dim);

      for i_func = 1:length(functions)
        func = functions(i_func);

        max_instance = min(size(results.mse{i_model, i_func, i_dim}, 1), find(~all(isnan(results.mse{i_model, i_func, i_dim}), 2), 1, 'last'));
        if (isempty(max_instance))
          continue;
        end

        for i_instance = 1:max_instance
          inst = instances(i_instance);

          for i_snapshot = 1:length(snapshots)
            snapshot = snapshots(i_snapshot);            
            allRow = allRow + 1;
            resultRowNamesAll{allRow} = sprintf('%s_%s_%dD_f%d_i%d_S%d', modelName, hash, dim, func, inst, snapshot);
            resultCellAll(allRow, 1:length(descrColNamesAll)) = {hash, dim, func, inst, snapshot};
          end
        end

        thisRows = lastFilledRow + (1:max_instance*length(snapshots));
        lastFilledRow = thisRows(end);
        for i_stat = 1:length(opts.statistics)
          stat = opts.statistics{i_stat};
          resultMatAll(thisRows, i_stat) = reshape(results.(stat){i_model, i_func, i_dim}, [], 1);
        end
        i_stat = length(opts.statistics)+1;
        if (~isempty(results.isTrained{i_model, i_func, i_dim}))
          resultMatAll(thisRows, i_stat) = reshape(results.isTrained{i_model, i_func, i_dim}, [], 1);
        else
          resultMatAll(thisRows, i_stat) = reshape(~isnan(results.mse{i_model, i_func, i_dim}), [], 1);
        end
        i_stat = i_stat + 1;
        if (isfield(results, 'isTrained2'))
          resultMatAll(thisRows, i_stat) = reshape(results.isTrained2{i_model, i_func, i_dim}, [], 1);
        end

        for agSnap = 1:length(opts.aggSnapshots)
          i_stat = 0;
          thisAggSnapshots = opts.aggSnapshots{agSnap};
          aggRow = aggRow + 1;
          resultRowNamesAgg{aggRow} = sprintf('%s_%s_%dD_f%d_a%d', modelName, hash, dim, func, agSnap);
          resultCellAgg(aggRow, 1:length(descrColNamesAgg)) = {hash, dim, func, agSnap};
          for is = 1:length(opts.statistics)
            i_stat = i_stat + 1;
            stat = opts.statistics{i_stat};
            resultMatAgg(aggRow, i_stat) = opts.aggFunction(results.(stat){i_model, i_func, i_dim}(:,thisAggSnapshots));
          end
          i_stat = i_stat + 1;
          if (~isempty(results.isTrained{i_model, i_func, i_dim}))
            resultMatAgg(aggRow, i_stat) = sum(sum(results.isTrained{i_model, i_func, i_dim}(:,thisAggSnapshots)));
          else
            resultMatAgg(aggRow, i_stat) = sum(sum(~isnan(results.mse{i_model, i_func, i_dim}(:,thisAggSnapshots))));
          end
          i_stat = i_stat + 1;
          if (isfield(results, 'isTrained2'))
            resultMatAgg(aggRow, i_stat) = sum(sum(results.isTrained2{i_model, i_func, i_dim}(:,thisAggSnapshots)));
          end
        end

      end  % for functions
    end  % for dimensions
  end  % for models

  % Final tables
  resultTableAll = cell2table(resultCellAll);
  resultTableAll(:, end+(1:size(resultMatAll,2))) = num2cell(resultMatAll);
  resultTableAll.Properties.VariableNames = resultColNamesAll;
  resultTableAll.Properties.RowNames = resultRowNamesAll;

  resultTableAgg = cell2table(resultCellAgg);
  resultTableAgg(:, end+(1:size(resultMatAgg,2))) = num2cell(resultMatAgg);
  resultTableAgg.Properties.VariableNames = resultColNamesAgg;
  resultTableAgg.Properties.RowNames = resultRowNamesAgg;

  writetable(resultTableAll, fullfile(exppath_short, exp_id, 'adaptModeltest_all_table.txt'), 'WriteRowNames', true);
  writetable(resultTableAgg, fullfile(exppath_short, exp_id, 'adaptModeltest_agg_table.txt'), 'WriteRowNames', true);

  save(fullfile(exppath_short, exp_id, ['modelStatistics_' datestr(now(), 'yyyy-mm-dd') '.mat']), 'resultTableAll', 'resultTableAgg', 'resultMatAll', 'resultMatAgg', 'folderModelOptions', 'modelFolders', 'functions', 'dimensions', 'snapshots', 'instances');
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
