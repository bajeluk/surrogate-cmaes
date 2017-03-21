function modelFolder = testModels(modelOptions, opts, funcToTest, dimsToTest, instToTest)
% testModels -- tests fitting power of models on DTS dataset.
%
% modelFolder = testModels(...) tests models and returns resulting names 
% of folders.
%
% Input:
%   modelOptions - model options to test | struct (cell-array of struct)
%   opts         - other options and settings
%       .modelType      - type of model (from ECFactory) | string
%       .exp_id         - exp_id of the _model testing_ experiment | string
%       .dataset        - filename (w/o extension) | string
%       .maxEvals       - maximal number of FE per dimensions to consider | integer
%       .exppath_short  - directory with experiments subdirectories | string
%       .statistics     - cell-array of statistics' names | cell array of strings
%       .rewrite_results - whether to rewrite already saved results | bool
%   funcToTest         - functions to test | array of integers
%   dimsToTest         - dimensions to test | array of integers
%   instToTest         - dimensions to test | array of integers
%
% Output:
%   modelFolder  - list of folders containing results | cell-array of 
%                  string
%
% See Also:
%   modelTestSet, testOneModel

  % Default input parameters settings
  if nargin < 5
    instToTest = [1];
    if nargin < 4
      dimsToTest = [2];
      if nargin < 3
        funcToTest = [1];
      end
    end
  end

  % Ensure input parameters as cell arrays
  modelType = opts.modelType;
  if ~iscell(modelType)         modelType = {modelType}; end
  if ~iscell(modelOptions)      modelOptions = {modelOptions}; end
  nModel = length(modelOptions);
  if length(modelType) == 1 && nModel > 1
    modelType = modelType(ones(nModel, 1));
  end
  assert(isnumeric(funcToTest), '''funcToTest'' has to be integer')
  assert(isnumeric(dimsToTest), '''dimsToTest'' has to be integer')
  assert(isnumeric(instToTest), '''instToTest'' has to be integer')

  % Default options in 'opts'
  opts.maxEvals = defopts(opts, 'maxEvals', 250);
  opts.exppath_short = defopts(opts, 'exppath_short', fullfile('exp', 'experiments'));
  opts.statistics = defopts(opts, 'statistics', { 'mse' });

  % Load the dataset and its parameters
  assert(isfield(opts, 'dataset'), 'opts does not include a dataset');
  if ischar(opts.dataset)
    assert(exist(opts.dataset, 'file')==2, 'Dataset %s does not exist', opts.dataset)
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
  maxEvals      = min(dataMaxEvals, opts.maxEvals);

  % filter the required functions/dimensions/instances to those that are really
  % in the dataset
  dimsToTest = restricToDataset(dimsToTest, dataDims, 'Dimensions');
  funcToTest = restricToDataset(funcToTest, dataFunc, 'Functions');
  instToTest = restricToDataset(instToTest, dataInst, 'Instances');

  % Assign (hopefully) unique names and output directories to models
  for m = 1:nModel
    modelHashName{m} = [modelType{m}, 'model_', modelHash(modelOptions{m})];
    % TODO: ensure, that this name is really unique
    modelFolder{m}   = fullfile(opts.exppath_short, opts.exp_id, [modelHashName{m}, '_', num2str(maxEvals), 'FE']);
  end

  % dimension loop
  for dim = dimsToTest
    d_data = find(dim == dataDims, 1);

    % function loop
    for fun = funcToTest
      f_data = find(fun == dataFunc, 1);

      % Check that we really have this data
      if (isempty(data{f_data, d_data}))
        warning(sprintf('Data of func %d in dim %d is missing in the dataset!', dataFunc(f_data), dataDims(d_data)));
        continue;
      end

      % model loop
      for m = 1:nModel
        fprintf('*******************  Fun: %d  Dim: %d  Model: %d  *******************\n', fun, dim, m)

        % directory and filename for the output
        [~, ~] = mkdir(modelFolder{m});
        modelFile = fullfile(modelFolder{m}, sprintf('%s_f%d_%dD.mat', modelHashName{m}, fun, dim));

        % do not rewrite existing files unless wanted
        if (exist(modelFile, 'file') && ~opts.rewrite_results)
          % there are some already calculated results, try to load them and continue
          oldResults = load(modelFile, 'stats', 'models', 'y_models', 'instances', 'modelOptions', 'fun', 'dim');
          finishedInstances = ismember(instToTest, oldResults.instances);
          if (all(finishedInstances))
            fprintf('All instances in %s already calculated. Skipping model testing.\n', modelFile);
            continue;
          else
            startInstanceIdx = find(~finishedInstances, 1);
            fprintf('Some instances already calculated in %s.\n', modelFile);
            fprintf('Starting from inst. # %d.\n', instToTest(startInstanceIdx));
            stats = oldResults.stats;
            models = oldResults.models;
            y_models = oldResults.y_models;
          end
        else
          % (re-)calculate the results
          startInstanceIdx = 1;
          if (exist(modelFile, 'file'))
            warning('Stop testing if you do not want to rewrite file %s', modelFile)
          end

          % prepare output variables
          stats = struct();
          for st = 1:length(opts.statistics)
            stats.(opts.statistics{st}) = NaN(length(instToTest), dataNSnapshots);
          end
          models   = cell(length(instToTest), dataNSnapshots);
          y_models = cell(length(instToTest), dataNSnapshots);
        end

        % instances loop
        for ii = startInstanceIdx:length(instToTest)
          inst = instToTest(ii);
          i_data = find(inst == dataInst, 1);
          fprintf('-- instance %2d --\n', inst);

          % train & test the model on the 'dataNSnapshots' datasets
          % from the current instance
          [new_stats, models(ii, :), y_models(ii, :)] = ...
              testOneModel(modelType{m}, modelOptions{m}, ...
              data{f_data, d_data, i_data}, dataNSnapshots, opts);

          % save results into output variables
          for st = 1:length(opts.statistics)
            fname = opts.statistics{st};
            stats.(fname)(ii, :) = new_stats.(fname);
          end

          % save results of the so-far calculated instances
          instances = instToTest(1:ii);
          save(modelFile, 'stats', 'models', 'y_models', 'instances', 'modelOptions', 'fun', 'dim')
        end  % instance loop
      end  % model loop
    end  % function loop
  end  % dimension loop

end

function hash = modelHash(modelOptions)
%TODO: proper model hash
% function creating hash for model identification using modelOptions
    if (isempty(modelOptions) || ~isstruct(modelOptions))
      hash = '0';
      return;
    end

    % gain fields and values of modelOptions
    [modelField, modelValues] = getFieldsVals(modelOptions);

    S = printStructure(modelOptions, 'Format', 'field');
    S = double(S);
    % exclude not necessary characters
    S = S(S > 32 & S~= 61) - 32;

    % create hash
    hash = num2str(sum(S.*(1:length(S))));
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
