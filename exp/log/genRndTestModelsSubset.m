function genRndTestModelsSubset(filename, varargin)
%genRndTestModelsSubset generate random subset of dataset devoted to model
% testing by testModels.
%
% genRndTestModelsSubset(filename, settings) random subset of dataset
% 'filename' will be generated according to 'settings'.
%
% Input:
%   filename - name of original dataset | string
%   settings  - pairs of property (string) and value, or struct with 
%               properties as fields:
%     'NSets'      - number of subsets to select | integer
%     'ResultName' - name of resulting file | string
%     'VarRange'   - ranges of variables (ie. which values to incorporate);
%                    empty for all values | 5x1 cell-array of double
%     'VarState'   - states of variables | 5x1 cell-array
%                      'rnd'  - randomly select values from 'VarRange'
%                      'full' - use all values from 'VarRange'
%
% See Also:
%   datasetFromInstances, modelTestSets, testModels

  if nargin < 1
    help getRndTestModelsSubset
    return
  end
  
  % check input
  assert(isfile(filename), 'There is no dataset file called %s', filename)
  
  % parse settings
  settings = settings2struct(varargin);
  [pathstr, name, ext] = fileparts(filename);
  resName = defopts(settings, 'ResultName', fullfile(pathstr, [name, '_subset', ext]));
  varRange = defopts(settings, 'VarRange', {});
  varState = defopts(settings, 'VarState', {'rnd', 'rnd', 'rnd', 'rnd', 'rnd'});
  nSets = defopts(settings, 'NSets', 100);
  
  % load file
  S = load(filename);
  
  % checkout ranges 
  % (be careful: generation range is in ids, not real numbers)
  defRanges = {S.fun, S.dim, S.inst, S.modelSettings, 1:S.nSnapshotsPerRun};
  if isempty(varRange)
    % all defaults
    varRange = defRanges;
  elseif any(cellfun(@isempty, varRange))
    % some defaults
    emptyId = cellfun(@isempty, varRange);
    varRange(emptyId) = defRanges(emptyId);
  end
  % TODO: check values if do they exist
  
  % generate ids
  isRnd = cellfun(@(x) strcmp(x, 'rnd'), varState);
  fullVec = varRange(~isRnd);
  % combine full variables
  fullId = combvec(fullVec{:})';
  % how many times to use each fullId row?
%   maxUsage = ceil(nSets/size(fullId, 1));
  rowIdUsage = cvInd(nSets, size(fullId, 1));
  % create final ids and select datasets
  ds = cell(1, nSets);
  vals = NaN(nSets, 5);
  for n = 1:nSets
    % add full variables
    vals(n, ~isRnd) = fullId(rowIdUsage(n), :);
    % random selection can be repeated if dataset does not exist
    selected = false;
    nSelect = 0;
    while ~selected && nSelect < numel(S.ds)*S.nSnapshotsPerRun
      nSelect = nSelect + 1;
      try
        % add random variables
        vals(n, isRnd) = cellfun(@(x) x(randi(numel(x))), varRange(isRnd));
        % get ids (except generations)
        ids = arrayfun(@(x) find(vals(n, x) == defRanges{x}, 1, 'first'), 1:4);
        % select appropriate dataset
        ds{n} = selectGen(S.ds{ids(1), ids(2), ids(3), ids(4)}, vals(n, 5));
        selected = true;
      catch
      end
    end
    % selection not successful
    if isempty(ds{n})
      warning('Dataset n. %d is empty', n)
    end
  end
  
  % extract unique values to save
  fun = unique(vals(:, 1));
  dim = unique(vals(:, 2));
  inst = unique(vals(:, 3));
  modelSettings = unique(vals(:, 4));
  generations = unique(vals(:, 5));
  
  % save results
  save(resName, 'ds', 'fun', 'dim', 'inst', 'modelSettings', 'generationIds')

end

function res = selectGen(dataset, genId)
% selectGen Selects generation from dataset

  fnames = fieldnames(dataset);
  for fn = 1:numel(fnames)
    actName = fnames{fn};
    if any(size(dataset.(actName)) > 1)
      % select appropriate generation according to its id
      if iscell(dataset.(actName))
        res.(actName) = dataset.(actName){genId};
      else
        res.(actName) = dataset.(actName)(genId);
      end
    else
      % copy field as it is
      res.(actName) = dataset.(actName);
    end
  end
end