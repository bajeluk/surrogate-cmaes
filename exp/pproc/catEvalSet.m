function [evals, settings] = catEvalSet(folders, funcSet)
% Finds unique settings and prepares its data for further processing.
% [evals, settings] = catEvalSet(folders, funcSet) returns cell array 
% 'data' of size functions x dimensions x unique settings and appropriate 
% 'settings'.
%
% Input:
%   folders - path to data | string or cell-array of strings
%   funcSet - structure with fields 'BBfunc' (numbers of BBOB functions) 
%             and 'dims' (numbers of dimensions) | structure
%
% Output:
%   evals    - aggregated data of size functions x dimensions x settings 
%              | cell array
%   settings - appropriate settings to 'data' | structure
%
% See Also:
%   catEvalSet
  
  evals = {};
  settings = {};
  if nargin < 2
    help catEvalSet
    return
  end
  % cell checking
  if ~iscell(folders)
    folders = {folders};
  end
  
  % initialize
  nFolders = length(folders);
  exp_evals = cell(1, nFolders);
  settings = cell(1, nFolders);

  % load data from all folders
  for s = 1:length(folders)
    [exp_evals{s}, settings{s}] = dataReady(folders{s}, funcSet);
  end
  
  % not empty settings
  notEmptySet = ~cellfun(@isempty, settings);
  % add extra settings to empty ones
  emptyID = inverseIndex(~notEmptySet);
  for i = emptyID
    [~, algName] = fileparts(folders{i});
    settings{i} = {struct('algName', algName)};
  end

  % gain all not empty settings
  allSettings = [settings{:}];
  % remove field 'experimentPath' because it is different for each
  % experiment
  expPathFieldID = cellfun(@(x) isfield(x, 'experimentPath'), allSettings);
  allSettings(expPathFieldID) = cellfun(@(x) rmfield(x, 'experimentPath'), allSettings(expPathFieldID), 'UniformOutput', false);
  % find unique settings
  %TODO: efective finding of unique settings and ID's. Sth like:
  % help_settings = settings;
  % notEmptySet = ~cellfun(@isempty, help_settings);
  % while any(notEmptySet)
  %   settingsID(getStructIndex(settings, settings{find(notEmptySet, 1, 'first')})) = s;
  %   help_settings(settingsID == s) = {};
  %   notEmptySet = ~cellfun(@isempty, help_settings);
  % end
  for s = length(allSettings):-1:1
    settingsID(getStructIndex(allSettings, allSettings{s})) = s;
  end
  
  % concatenate evaluations from different experiments and with the same
  % settings
  exp_evals = cat(3, exp_evals{:});
  nSettings = length(allSettings);
  evals = cell(length(funcSet.BBfunc), length(funcSet.dims), nSettings);
  for s = 1 : nSettings
    for f = 1 : length(funcSet.BBfunc)
      for d = 1 : length(funcSet.dims)
          evals{f, d, s} = [exp_evals{f, d, settingsID == s}];
      end
    end
  end
  
  % return unique settings and its evaluations
  settings = allSettings(unique(settingsID));
  evals = evals(:, :, unique(settingsID));

end