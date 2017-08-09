function [algData, algNames, algColors] = loadCompAlg(path, funcSet)
% [algData, algNames, algColors] = loadCompAlg(path, funcSet) loads data,
% names and colors of compared algorithms for functions and dimensions 
% given.

  if nargin < 1
    path = fullfile('exp', 'pproc', 'compAlgMat.mat');
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
  
  % extract required functions and dimensions
  nFuns = length(funcSet.BBfunc);
  nDims = length(funcSet.dims);
  for a = 1:length(algorithms)
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
    
  % gain names and colors
  algNames = {algorithms.name};
  algColors = cell2mat({algorithms.color}');

end

function updateAlg(path)
% updates file with algorithms' data

  % download latest version
  try
    websave(path, 'http://artax.karlin.mff.cuni.cz/~bajel3am/scmaes/compAlgMat.mat');
  catch
  end
end