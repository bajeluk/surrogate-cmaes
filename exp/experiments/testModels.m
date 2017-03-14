function modelFolder = testModels(modelType, modelOptions, opts, func, dims, instances)
% testModels -- tests fitting power of models on DTS dataset.
%
% modelFolder = testModels(...) tests models and returns resulting names 
% of folders.
%
% Input:
%   modelType    - type of model (from ECFactory) | string
%   modelOptions - model options to test | struct (cell-array of struct)
%   dataset      - data for testing | cell-array or string
%   func         - functions to test | double
%   dims         - dimensions to test | double
%   rewrite      - rewrite existing result files | boolean
%
% Output:
%   modelFolder  - list of folders containing results | cell-array of 
%                  string
%
% See Also:
%   modelTest, modelTestSet, Model

  % path settings
  defFolder = fullfile('exp', 'experiments', 'model');
  defDatasetPath = fullfile(defFolder, 'defData', 'defSet.mat');
 
  if nargin < 6
    if nargin < 5
      if nargin < 4
        if nargin < 3
          if nargin < 2
            if nargin < 1
              help testModels
              return
            end
            modelOptions = struct();
          end
          dataset = defDatasetPath;
        end
        func = 1;
      end
      dims = 2;
    end
    rewrite = false;
  end
  
  % input check
  if ~iscell(modelType)
    modelType = {modelType};
  end
  if ~iscell(modelOptions)
    modelOptions = {modelOptions};
  end
  nModel = length(modelOptions);
  if length(modelType) == 1 && nModel > 1
    modelType = modelType(ones(nModel, 1));
  end
  if ischar(dataset)
    assert(exist(dataset, 'file')==2, 'Dataset %s does not exist', dataset)
  else
    assert(iscell(dataset), 'Dataset has to be string or cell')
  end
  assert(isnumeric(func), '''func'' has to be integer')
  assert(isnumeric(dims), '''dims'' has to be integer')
  
  % load data
  if ischar(dataset)
    loadData = load(dataset);
  else
    loadData = struct();
  end
  data = defopts(loadData, 'ds', dataset);
  dataFun = defopts(loadData, 'fun', func);
  dataDim = defopts(loadData, 'dim', dims);
  maxEval = defopts(loadData, 'maxEval', []);
  
  % check functions and dimensions for testing
  funcId = ismember(func, dataFun);
  dimId = ismember(dims, dataDim);
  if ~all(funcId)
    warning('Functions %s are not in tested dataset', ...
      strjoin(arrayfun(@num2str, func(~funcId), 'UniformOutput', false), ','))
    func = func(funcId);
  end
  if ~all(dimId)
    warning('Dimensions %s are not in tested dataset', ...
      strjoin(arrayfun(@num2str, dims(~dimId), 'UniformOutput', false), ','))
    dims = dims(dimId);
  end
  nFunc = length(func);
  nDims = length(dims);
  
  modelName = cellfun(@(x,y) [x, 'Model_', modelHash(y)], modelType, modelOptions, 'UniformOutput', false);
  modelFolder = cellfun(@(x) fullfile(defFolder, x), modelName, 'UniformOutput', false);
    % if there is an information about maximum of evaluations write it
  if ~isempty(maxEval)
    modelFolder = cellfun(@(x) [x, '_', num2str(maxEval), 'FE'], modelFolder, 'UniformOutput', false);
  end
  for m = 1:nModel    
    if ~exist(modelFolder{m}, 'dir')
      mkdir(modelFolder{m})
    end
  end
  
  % dimension loop
  for dim = dims
    d_data = find(dim == dataDim, 1);
    % function loop
    for fun = func
      f_data = find(fun == dataFun, 1);
      if ~isempty(data{f_data, d_data})
        % model loop
        for m = 1:nModel
          fprintf('*******************  Fun: %d  Dim: %d  Model: %d  *******************\n', fun, dim, m)

          modelFile = fullfile(modelFolder{m}, sprintf('%s_f%d_%dD.mat', modelName{m}, fun, dim));
          % do not rewrite existing files unless wanted
          if exist(modelFile, 'file') && ~rewrite
            fprintf('File %s already exist.\n', modelFile)
          else
            % warn user if the result file will replaced by the new one
            if exist(modelFile, 'file') && rewrite
              warning('Stop testing if you do not want to rewrite file %s', modelFile)
            end
            % test model
            [stats, model, ym] = modelTest(modelType{m}, modelOptions{m}, data{f_data, d_data});
            % save model results
            save(modelFile, 'stats', 'model', 'ym')
          end
        end
      end
    end
  end
  
end

function hash = modelHash(modelOptions)
%TODO: proper model hash
% function creating hash for model identification using modelOptions

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
