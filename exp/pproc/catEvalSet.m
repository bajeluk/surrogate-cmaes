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
  % remove empty settings and appropriate exp_evals
  notEmptySet = ~cellfun(@isempty, settings);
  if ~any(notEmptySet)
    evals = exp_evals;
    return
  end
  settings = settings(notEmptySet);
  exp_evals = exp_evals(notEmptySet);
  % remove field 'experimentPath' because it is different for each
  % experiment
  allSettings = [settings{:}];
  expPathFieldID = cellfun(@(x) isfield(x, 'experimentPath'), allSettings);
  settings = cellfun(@(x) rmfield(x, 'experimentPath'), allSettings(expPathFieldID), 'UniformOutput', false);
  % find unique settings
  %TODO: efective finding of unique settings and ID's. Sth like:
  % help_settings = settings;
  % notEmptySet = ~cellfun(@isempty, help_settings);
  % while any(notEmptySet)
  %   settingsID(getStructIndex(settings, settings{find(notEmptySet, 1, 'first')})) = s;
  %   help_settings(settingsID == s) = {};
  %   notEmptySet = ~cellfun(@isempty, help_settings);
  % end
  for s = length(settings):-1:1
    settingsID(getStructIndex(settings, settings{s})) = s;
  end
  settings = settings(unique(settingsID));
  
  % concatenate evaluations from different experiments and with the same
  % settings
  exp_evals = cat(3, exp_evals{:});
  nSettings = length(settings);
  evals = cell(length(funcSet.BBfunc), length(funcSet.dims), nSettings);
  for s = 1 : nSettings
    for f = 1 : length(funcSet.BBfunc)
      for d = 1 : length(funcSet.dims)
        evals{f, d, s} = [exp_evals{f, d, settingsID == s}];
      end
    end
  end

end