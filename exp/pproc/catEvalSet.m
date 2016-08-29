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
  notEmptyID = find(notEmptySet);
  settings_ne = settings(notEmptySet);
  % not empty evaluations
  exp_evals_ne = exp_evals(notEmptySet);
  % remove field 'experimentPath' because it is different for each
  % experiment
  allSettings = [settings_ne{:}];
  expPathFieldID = cellfun(@(x) isfield(x, 'experimentPath'), allSettings);
  settings_ne = cellfun(@(x) rmfield(x, 'experimentPath'), allSettings(expPathFieldID), 'UniformOutput', false);
  % find unique settings
  %TODO: efective finding of unique settings and ID's. Sth like:
  % help_settings = settings;
  % notEmptySet = ~cellfun(@isempty, help_settings);
  % while any(notEmptySet)
  %   settingsID(getStructIndex(settings, settings{find(notEmptySet, 1, 'first')})) = s;
  %   help_settings(settingsID == s) = {};
  %   notEmptySet = ~cellfun(@isempty, help_settings);
  % end
  for s = length(settings_ne):-1:1
    settingsID(getStructIndex(settings_ne, settings_ne{s})) = s;
  end
  
  % concatenate evaluations from different experiments and with the same
  % settings
  exp_evals_ne = cat(3, exp_evals_ne{:});
  nSettings = length(settings_ne);
  evals = cell(length(funcSet.BBfunc), length(funcSet.dims), nSettings);
  for s = 1 : nSettings
    for f = 1 : length(funcSet.BBfunc)
      for d = 1 : length(funcSet.dims)
        if ~isempty(settings_ne{s})
          evals{f, d, s} = [exp_evals_ne{f, d, settingsID == s}];
        else
          evals{f, d, s} = exp_evals{s}{f, d};
        end
      end
    end
  end
  
  % return unique settings and its evaluations
  % [folders with settings, folders without settings]
  settings = [settings_ne(unique(settingsID)), folders(~notEmptySet)];
  evals = cat(3, evals(:,:,ismember(1:nSettings, unique(settingsID))), evals(:, :, ~notEmptySet));

end