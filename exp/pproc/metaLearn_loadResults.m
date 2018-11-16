function [results, settings, resParams] = metaLearn_loadResults(folders, varargin)
% [results, settings, resParams] = metaLearn_loadResults(folders, varargin) 
% load results of metalearning experiments.
%
% Input:
%   folders  - path to data | string or cell-array of strings
%   varargin - pairs of property (string) and value or struct with 
%              properties as fields:
%     'IgnoreSettings' - list of fields in meta-settings to ignore (e.g. in
%                        case of random initialization) | string or 
%                        cell-array of strings
%     'ShowOutput'     - print output of data loading to screen | boolean
%
% Output:
%   results   - aggregated results of size functions x dimensions x 
%               settings | cell-array of structures
%   settings  - appropriate settings to results | cell-array of struct
%   resParams - parameters of loaded results (and experiment settings) |
%               struct
%
% See Also:
%   dataReady, catEvalSet

  if nargout > 0
    results = {};
    settings = {};
  end
  if nargin < 1
    help metaLearn_loadResults
    return
  end
  % cell checking
  if ~iscell(folders)
    folders = {folders};
  end
  function_settings = settings2struct(varargin);
  
  % parse settings
  

  % load data from first folder (temporary)
  [results, settings, resParams] = metaLearn_loadOneFolder(folders{1}, function_settings);
  
  % load data from all folders
%   for s = 1:length(folders)
%     [exp_results{s}, settings{s}] = metaLearn_loadOneFolder(folders{s});
%   end
end

function [results, settings, res_params] = metaLearn_loadOneFolder(exp_folder, fSettings)
% load results from one folder
  
  ignoreSettings = defopts(fSettings, 'IgnoreSettings', {''});
  if ~iscell(ignoreSettings)
    ignoreSettings = {ignoreSettings};
  end
  showOutput = defopts(fSettings, 'ShowOutput', false);

  % load parameter file
  paramFile = fullfile(exp_folder, 'metalearn_params.mat');
  if isfile(paramFile)
    params = load(paramFile);
    if showOutput
      fprintf('Experiment setting file loaded\n')
    end
  else
    warning('Experiment setting file is missing (metalearn_params.mat)')
    params = [];
  end
  % gain file list
  if showOutput
    fprintf('Searching for MAT-files in folder %s\n', exp_folder)
  end
  fileList = searchFile(exp_folder, '*.mat');
  % remove param file
  if ~isempty(params)
    paramsId = strcmp(paramFile, fileList);
    fileList(paramsId) = [];
  end
  
  settings = {};
  results = {};
  instances = {};
  functions = [];
  dimensions = [];
  nFiles = numel(fileList);
  % load file sequentially
  for f = 1:nFiles
    if showOutput
      fprintf('Loading %s\n', fileList{f})
    else
      warning('off', 'MATLAB:load:variableNotFound')
    end
    S = load(fileList{f}, '-mat', 'dim', 'func', 'inst_samp', 'modelOpt', 'modelType', 'results', 'cv');
    warning('on', 'MATLAB:load:variableNotFound')
    if all(isfield(S, {'dim', 'func', 'inst_samp', 'modelOpt', 'modelType', 'results', 'cv'}))
      % unify parameters to one settings structure
      actualSettings = S.modelOpt;
      actualSettings.modelType = S.modelType;
      % if isfield(S, 'inst_samp')
      %   actualSettings.inst_samp = S.inst_samp;
      % end
      
      % ignore chosen settings
      if ~isempty(ignoreSettings{1})
        actualSettings = safermfield(actualSettings, ignoreSettings);
      end
      
      % load dimension, function, and settings
      [dimId, dimensions] = propertyId(dimensions, S.dim);
      [funcId, functions] = propertyId(functions, S.func);
      [settingsId, settings] = propertyId(settings, actualSettings);
      
      % save results to proper position
      results{funcId, dimId, settingsId} = S.results;
      % save instances to proper position
      instances{funcId, dimId, settingsId} = S.inst_samp;

    end
  end
  
  % sort according to function and dimension 
  [functions, funSortId] = sort(functions);
  [dimensions, dimSortId] = sort(dimensions);
  results = results(funSortId, dimSortId, :);
  instances = instances(funSortId, dimSortId, :);
  
  % return loaded functions, dimenstion and experiment settings
  res_params.functions = functions;
  res_params.dimensions = dimensions;
  res_params.metalearn_params = params;
  res_params.instances = instances;
  
end

function [propId, property] = propertyId(property, value)
% Gain property identifier. If the property value is new, add it to
% property list.
  if iscell(property)
    propId = cellfun(@(x) isequal(value, x), property);
    if isempty(propId) || ~any(propId)
      property{end+1} = value;
      propId = length(property);
    else
      propId = find(propId);
    end
  else
    propId = property == value;
    if isempty(propId) || ~any(propId)
      property(end+1) = value;
      propId = length(property);
    else
      propId = find(propId);
    end
  end
end