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
  % find equal settings
  nAllSettings = length(allSettings);
  equalSettings = logical(eye(nAllSettings));
  for s = 1 : (nAllSettings - 1)
    for t = (s + 1) : nAllSettings
      equalSettings(s, t) = isequal(allSettings{s}, allSettings{t});
    end
  end
  % find unique settings
  notUsed = true(1, nAllSettings);
  s = 0;
  while any(notUsed)
    s = s + 1;
    r = find(notUsed, 1, 'first');
    settingsID(notUsed & equalSettings(r, :)) = s;
    notUsed = notUsed & ~equalSettings(r, :);
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
  outputSettingsId = arrayfun(@(x) find(settingsID == x, 1, 'first'), unique(settingsID));
  settings = allSettings(outputSettingsId);
  evals = evals(:, :, unique(settingsID));

end