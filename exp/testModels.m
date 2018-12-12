function modelFolder = testModels(modelOptions, opts, funcToTest, dimsToTest, instToTest, idsToTest)
% testModels -- tests fitting power of models on DTS dataset.
%
% modelFolder = testModels(...) tests models and returns resulting names
% of folders.
%
% Input:
%   modelOptions - model options to test | struct (cell-array of struct)
%   opts         - other options and settings | struct
%     Compulsory fields:
%       .dataset         - filename (w/o extension) | string
%       .exp_id          - exp_id of the _model testing_ experiment |
%                          string
%       .modelType       - type of model (from ECFactory) | string
%     Optional fields:
%       .alwaysRetrain   - let always be the models be retrained
%                          (efectively deletes models in the 'data')
%       .exppath_short   - directory with experiments subdirectories |
%                          string
%       .maxEvals        - maximal number of FE per dimensions to consider
%                          | integer
%       .mfts_settings   - settings of data metafeatures (keep empty if no
%                          metafeature calculation should be performed) |
%                          struct
%       .rewrite_results - whether to rewrite already saved results | bool
%       .saveModels      - whether the models should be saved, too, default
%                          FALSE | bool
%       .scratch         - scrath directory | string
%       .statistics      - cell-array of statistics' names | cell array of
%                          strings
%   funcToTest   - functions  to test | array of integers
%   dimsToTest   - dimensions to test | array of integers
%   instToTest   - instances  to test | array of integers
%   idsToTest    - ids to test | array of integers
%
% Output:
%   modelFolder  - list of folders containing results | cell-array of
%                  string
%
% See Also:
%   modelTestSets, testOneModel, getDataMetaFeatures

  % Default input parameters settings
  if nargin < 6
    idsToTest = 1;
    if nargin < 5
      instToTest = 1;
      if nargin < 4
        dimsToTest = 2;
        if nargin < 3
          funcToTest = 1;
          if nargin < 2
            if nargout > 0
              modelFolder = {};
            end
            help testModels
            return
          end
        end
      end
    end
  end

  % Ensure input parameters as cell arrays
  assert(isfield(opts, 'modelType'), 'opts does not include a modelType field');
  assert(isfield(opts, 'exp_id'), 'opts does not include a experiment ID (exp_id)');
  modelType = opts.modelType;
  if ~iscell(modelType),        modelType = {modelType}; end
  if ~iscell(modelOptions),     modelOptions = {modelOptions}; end
  nModel = length(modelOptions);
  if length(modelType) == 1 && nModel > 1
    modelType = modelType(ones(nModel, 1));
  end
  assert(isnumeric(funcToTest), '''funcToTest'' has to be integer')
  assert(isnumeric(dimsToTest), '''dimsToTest'' has to be integer')
  assert(isnumeric(instToTest), '''instToTest'' has to be integer')
  assert(isnumeric(idsToTest), '''idsToTest'' has to be integer')

  % Default options in 'opts'
  opts.maxEvals = defopts(opts, 'maxEvals', 250);
  opts.exppath_short = defopts(opts, 'exppath_short', fullfile('exp', 'experiments'));
  opts.statistics = defopts(opts, 'statistics', { 'mse' });
  opts.saveModels = defopts(opts, 'saveModels', false);
  opts.scratch    = defopts(opts, 'scratch', fullfile(opts.exppath_short, opts.exp_id));
  opts.rewrite_results = defopts(opts, 'rewrite_results', false);
  opts.mfts_settings = defopts(opts, 'mfts_settings', []);

  % Load the dataset and its parameters
  assert(isfield(opts, 'dataset'), 'opts does not include a dataset');
  if ischar(opts.dataset)
    assert(exist(opts.dataset, 'file')==2, 'Dataset %s does not exist', opts.dataset)
    fprintf('Loading %s\n', opts.dataset)
    loadData = load(opts.dataset);
  else
    assert(iscell(opts.dataset), 'Dataset has to be string or cell')
    loadData = struct();
  end
  data          = defopts(loadData, 'ds', opts.dataset);
  dataFunc      = defopts(loadData, 'fun', funcToTest);
  dataDims      = defopts(loadData, 'dim', dimsToTest);
  dataInst      = defopts(loadData, 'inst', instToTest);
  dataMaxEvals  = defopts(loadData, 'maxEval', 250);
  dataNSnapshots = defopts(loadData, 'nSnapshotsPerRun', 10);
  dataIds       = defopts(loadData, 'modelSettings', idsToTest);
  maxEvals      = min(dataMaxEvals, opts.maxEvals);

  % filter the required functions/dimensions/instances to those that are really
  % in the dataset
  dimsToTest = restricToDataset(dimsToTest, dataDims, 'Dimensions');
  funcToTest = restricToDataset(funcToTest, dataFunc, 'Functions');
  instToTest = restricToDataset(instToTest, dataInst, 'Instances');
  idsToTest = restricToDataset(idsToTest, dataIds, 'Model settings');

  % Assign (hopefully) unique names and output directories to models
  for m = 1:nModel
    modelHashName{m} = [modelType{m}, 'model_', modelHash(modelOptions{m})];
    % TODO: ensure, that this name is really unique
    modelFolder{m}   = fullfile(opts.exppath_short, opts.exp_id, [modelHashName{m}, '_', num2str(maxEvals), 'FE']);
  end
  % create folder for metafeatures
  mftsFolder = fullfile(opts.exppath_short, opts.exp_id, 'metafeatures');
  [~, ~] = mkdir(mftsFolder);

  % dimension loop
  for dim = dimsToTest
    d_data = find(dim == dataDims, 1);

    % function loop
    for fun = funcToTest
      f_data = find(fun == dataFunc, 1);

      % Check that we really have this data
      if (isempty(data{f_data, d_data}))
        warning('Data of func %d in dim %d is missing in the dataset!', dataFunc(f_data), dataDims(d_data));
        continue;
      end

      % model loop
      for m = 1:nModel
        fprintf('*******************  Fun: %d  Dim: %d  Model: %d  *******************\n', fun, dim, m)

        % directory and filename for the output
        [~, ~] = mkdir(modelFolder{m});
        modelFile = fullfile(modelFolder{m}, sprintf('%s_f%d_%dD.mat', modelHashName{m}, fun, dim));
        missingDataFile = fullfile(modelFolder{m}, sprintf('%s_f%d_%dD_missing_data.mat', modelHashName{m}, fun, dim));

        % do not rewrite existing files unless wanted
        finishedInstances = false(size(instToTest));
        finishedIds = false(size(idsToTest));
        if (exist(modelFile, 'file') && ~opts.rewrite_results)
          % there are some already calculated results, try to load them and continue
          try
            if (opts.saveModels)
              oldResults = load(modelFile, 'stats', 'models', 'y_models', 'instances', 'ids', 'modelOptions', 'fun', 'dim');
            else
              oldResults = load(modelFile, 'stats', 'y_models', 'instances', 'ids', 'modelOptions', 'fun', 'dim');
            end
            finishedInstances = ismember(instToTest, oldResults.instances);
            finishedIds = ismember(idsToTest, oldResults.ids);
          catch err
            warning('File %s is corrupted. Will be rewritten. Error: %s', modelFile, err.message);
          end
          if (all(finishedInstances) && all(finishedIds))
            fprintf('All instances in %s already calculated. \nSkipping model testing.\n', modelFile);
            continue;
          end
        end

        % matrix of missing dim, fun, inst, id combinations
        if (exist(missingDataFile))
          load(missingDataFile, 'missingData');
        else
          missingData = zeros(0, 4);
        end

        % prepare output variables
        stats = struct();
        for st = 1:length(opts.statistics)
          stats.(opts.statistics{st}) = NaN(length(instToTest), length(idsToTest), dataNSnapshots);
        end
        models   = cell(length(instToTest), length(idsToTest), dataNSnapshots);
        y_models = cell(length(instToTest), length(idsToTest), dataNSnapshots);

        if (any(finishedInstances) && any(finishedIds))
          fprintf('Some instances and ids already calculated in %s.\n', modelFile);
          % copy in finished hypercube
          for i_inst = 1:length(instToTest)
            for i_id = 1:length(idsToTest)
              if finishedInstances(i_inst) && finishedIds(i_id)
                j_inst = find(instToTest(i_inst) == oldResults.instances, 1);
                j_id = find(idsToTest(i_id) == oldResults.ids, 1);
                for st = 1:length(opts.statistics)
                  stats.(opts.statistics{st})(i_inst, i_id, :) = oldResults.stats.(opts.statistics{st})(j_inst, j_id, :);
                end

                if (opts.saveModels)
                  models{i_inst, i_id} = oldResults.models{j_inst, j_id};
                end
                y_models{i_inst, i_id} = oldResults.y_models{j_inst, j_id};
              end
            end
          end
        else
          % (re-)calculate the results
          if (exist(modelFile, 'file'))
            warning('Stop testing if you do not want to rewrite file %s', modelFile)
          end
        end

        % print what you are going to test
        printStructure(modelOptions{m}, 'StructName', 'modelOptions')

        % instances loop
        for ii = 1:length(instToTest)
          for i_id = 1:length(idsToTest)
            inst = instToTest(ii);
            id = idsToTest(i_id);
            i_data = find(inst == dataInst, 1);
            id_data = find(id == dataIds, 1);

            idx = [dim, fun, inst, id];
            is_missing = ismember(missingData, idx, 'rows');
            if isempty(data{f_data, d_data, i_data, id_data})
              fprintf('Missing data %d %d %d %d\n', dim, fun, inst, id);
              if ~any(is_missing)
                missingData(end+1, :) = idx;
              end
              save(missingDataFile, 'missingData');

              continue;
            else
              % mark indices as not missing
              if any(is_missing)
                fprintf('Formerly missing data are not missing %d %d %d %d\n', dim, fun, inst, id);
                missingData(is_missing, :) = [];
                save(missingDataFile, 'missingData');
              end
            end

            if ~any(is_missing) && finishedInstances(ii) && finishedIds(i_id)
              fprintf('Skipping already finished inst %d, id %d\n', ii, i_id);
              continue;
            end

            fprintf('-- instance %2d, id %2d --\n', inst, id);

            BBOB_DIR = [opts.scratch '/bbob_output'];
            fgeneric('initialize', fun, inst, BBOB_DIR);
            % Info about tested BBOB function
            bbob_handlesF = benchmarks('handles');
            noisyHandles = benchmarksnoisy('handles');
            bbob_handlesF(100+(1:length(noisyHandles))) = noisyHandles;
            opts.bbob_func = bbob_handlesF{fun};

            % delete the saved models from the dataset to be always retrained
            if (isfield(opts, 'alwaysRetrain') && opts.alwaysRetrain)
              if (isfield(data{f_data, d_data, i_data, id_data}, 'models'))
                data{f_data, d_data, i_data, id_data}.models = cell(size(data{f_data, d_data, i_data, id_data}.models));
              end
              if (isfield(data{f_data, d_data, i_data, id_data}, 'models2'))
                data{f_data, d_data, i_data, id_data}.models2 = cell(size(data{f_data, d_data, i_data, id_data}.models2));
              end
            end

            % train & test the model on the 'dataNSnapshots' datasets
            % from the current instance
            if (~isfield(opts, 'trySecondModel') || ~opts.trySecondModel)
              [new_stats, thisModels, y_models(ii, i_id, :)] = ...
                  testOneModel(modelType{m}, modelOptions{m}, ...
                  data{f_data, d_data, i_data, id_data}, dataNSnapshots, opts);
              modelsVarString = { 'y_models' };
              if (opts.saveModels)
                modelsVarString(end+1) = 'models';
              end
            else
              [new_stats, thisModels, y_models(ii, i_id, :), thisModels2, y_models2(ii, i_id, :)] = ...
                  testOneModel(modelType{m}, modelOptions{m}, ...
                  data{f_data, d_data, i_data, id_data}, dataNSnapshots, opts);
              modelsVarString = { 'y_models', 'y_models2' };
              if (opts.saveModels)
                modelsVarString(end+(1:2)) = { 'models', 'models2' };
              end
            end

            % finalize BBOB function settings
            fgeneric('finalize');

            % save results into output variables
            for st = 1:length(opts.statistics)
              fname = opts.statistics{st};
              stats.(fname)(ii, i_id, :) = new_stats.(fname);
            end

            % save results of the so-far calculated instances
            instances = instToTest(1:ii);
            ids = idsToTest(1:i_id);
            if (opts.saveModels)
              models(ii, i_id, :) = thisModels;
              if (isfield(opts, 'trySecondModel') && opts.trySecondModel)
                models2(ii, i_id, :) = thisModels2;
              end
            end
            save(modelFile, 'stats', modelsVarString{:}, 'instances', 'ids', 'modelOptions', 'fun', 'dim');

            % calculate metafeatures of the dataset if not calculated earlier
            % or rewrite results setting is on
            opts.mfts_settings.fun = fun;
            opts.mfts_settings.dim = dim;
            opts.mfts_settings.instances = inst;
            opts.mfts_settings.output = ...
              sprintf('%s%sdata_f%d_%dD_inst%d_id%d_fts.mat', ...
                      mftsFolder, filesep, fun, dim, inst, id);
            if ~isempty(opts.mfts_settings) && ...
                (opts.rewrite_results || ~isfile(opts.mfts_settings.output))
              getDataMetaFeatures(data{f_data, d_data, i_data, i_id}, opts.mfts_settings)
            end
          end % id loop
        end  % instance loop
      end  % model loop
    end  % function loop
  end  % dimension loop

  fprintf('****************************  FINISH  ****************************\n')
end


function [sField, sVal] = getFieldsVals(s)
% sf = getFields(s, fields) extracts fields and its values from structure s

  sField = fieldnames(s);
  nFields = length(sField);
  sVal = cell(nFields, 1);
  for i = 1:nFields
    sVal{i} = s.(sField{i});
  end
end


function [intToTest] = restricToDataset(intToTest, dataInt, identifier)
  id_ok = ismember(intToTest, dataInt);
  if ~all(id_ok)
    warning('%s %s are not in tested dataset', identifier, ...
      strjoin(arrayfun(@num2str, intToTest(~id_ok), 'UniformOutput', false), ','))
    intToTest = intToTest(id_ok);
  end
end
