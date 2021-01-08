function [algData, algNames, algColors] = loadCompAlg(path, funcSet, algs)
% [algData, algNames, algColors] = loadCompAlg(path, funcSet, algs) loads
% data names and colors of compared algorithms for functions and dimensions
% given.
%
% Input:
%   path    - path to file containing reference algorithms' data | string |
%             default: 'exp/pproc/compAlgMat'
%   funcSet - structure defining functions and dimensions to load with
%             fields:
%               'BBfunc' - BBOB COCO function numbers | integer | default:
%                          1:24
%               'dims'   - function dimensions | integer | default: [2, 3,
%                          5, 10, 20]
%   algs    - names of algorithms to load | cell-array of following
%             strings:
%               'cmaes' - CMA-ES v. 3.61beta
%               'scmaes' - S-CMA-ES
%               'dtscmaes' - DTS-CMA-ES
%               'sacmes'   - BIPOP-s*ACMES-k (results from COCO web pages
%                            to paper loshchilov2013bipopulation)
%               'lmmcmaes' - lmm-CMA-ES (results from COCO web pages to
%                            paper auger2013benchmarking)
%               'smac'     - SMAC-BBOB (results from COCO web pages to
%                            paper hutter2013evaluation)
%           -  default: '' - loads all algorithms
%
% Output:
%   algData   - data for each algorithm | Nx1 cell-array
%   algNames  - names of individual loaded algorithms | Nx1 cell-array of
%               string
%   algColors - colors devoted to N loaded algorithms | Nx3 double

  if nargin < 3
    if nargin < 2
      if nargin < 1
        path = fullfile('exp', 'pproc', 'compAlgMat.mat');
      end
      funcSet.BBfunc = 1:24;
      funcSet.dims = [2, 3, 5, 10, 20];
    end
    algs = '';
  end
  algData = [];

  % update if does not exist
  if ~exist(path, 'file')
    updateAlg(path)
  end
  % update if cannot be loaded
  try
    alg = load(path);
  catch
    updateAlg(path)
    % load again
    try
      alg = load(path);
    % return if still cannot be loaded
    catch
      return
    end
  end
  
  % update if some fields are missing
  algorithms = alg.algorithm;
  if (~isfield(algorithms, 'BBfunc') || ~isfield(algorithms, 'dims'))
    updateAlg(path)
    % load again
    try
      alg = load(path);
      algorithms = alg.algorithm;
      assert(isfield(algorithms, 'BBfunc') || ~isfield(algorithms, 'dims'))
    % return if still cannot be loaded or fields are missing
    catch
      return
    end
  end
  
  % get proper algorithm names
  algs = properAlgNames(algs);
  algNames = {algorithms.name};
  algId = ismember(algNames, algs);

  % extract required functions and dimensions
  nFuns = length(funcSet.BBfunc);
  nDims = length(funcSet.dims);
  for a = find(algId)
    for d = 1:nDims
      dId = find(algorithms(a).dims == funcSet.dims(d));
      if ~isempty(dId)
        for f = 1:nFuns
          fId = find(algorithms(a).BBfunc == funcSet.BBfunc(f));
          if ~(isempty(fId) || isempty(dId))
            algData{a}{f,d} = algorithms(a).data{fId, dId};
          else
            algData{a}{f,d} = [];
          end
        end
      else
        algData{a}(1:nFuns,d) = cell(nFuns, 1);
      end
    end
  end
  % return only non-empty cells
  algData = algData(~cellfun(@isempty, algData));
    
  % gain names and colors
  algNames = algNames(algId);
  algColors = cell2mat({algorithms(algId).color}');

end

function updateAlg(path)
% updates file with algorithms' data

  % download latest version
  try
    websave(path, 'http://artax.karlin.mff.cuni.cz/~bajel3am/scmaes/compAlgMat.mat');
  catch
  end
end

function properNames = properAlgNames(codeNames)
% replace algorithm code names by proper ones

  defNames = {'CMA-ES', 'S-CMA-ES', 'DTS-CMA-ES', ...
              'BIPOP-{}^{s*}ACMES-k', 'SMAC', 'lmm-CMA-ES'};

  if ~iscell(codeNames)
    codeNames = {codeNames};
  end
  % if codenames are empty use all algorithms
  if numel(codeNames) == 1 && isempty(codeNames{1})
    properNames = defNames;
    return
  end

  properNames = cell(1, numel(codeNames));
  for a = 1:numel(codeNames)
    switch codeNames{a}
      case 'cmaes'
        properNames{a} = defNames{1};
      case 'scmaes'
        properNames{a} = defNames{2};
      case 'dtscmaes'
        properNames{a} = defNames{3};
      case 'sacmes'
        properNames{a} = defNames{4};
      case 'smac'
        properNames{a} = defNames{5};
      case 'lmmcmaes'
        properNames{a} = defNames{6};
      otherwise
        warning('Unknown algorithm codename ''%s''', codeNames{a})
        properNames{a} = '';
    end
  end
end