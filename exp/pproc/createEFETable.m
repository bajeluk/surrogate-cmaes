function [rankTable, ranks] = createEFETable(data, varargin)
% [table, ranks] = createEFETable(data, settings)
% Creates table containing rankings of expected number of function 
% evaluations (EFE) for chosen quantiles.
%
% Input:
%   data      - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataDims'  - dimensions of data
%     'DataFuns'  - functions of data
%     'Quantiles' - quantiles chosen to count
%     'TableDims' - dimensions chosen to count
%     'TableFuns' - functions chosen to count
%     'Target'    - EFE target
%
% Output:
%   table - table of EFE rankings
%   ranks - EFE rankings for each function and dimension
%
% See Also:
%   rankingTable, speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  rankTable = [];
  if nargin < 1 || isempty(data)
    help createEFETable
    return
  end
  settings = settings2struct(varargin{:});
  
  numOfData = length(data);
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  quantiles = defopts(settings, 'Quantiles', [0.25 0.5 0.75]);
  dims    = defopts(settings, 'TableDims', funcSet.dims);
  BBfunc  = defopts(settings, 'TableFuns', funcSet.BBfunc);
  efeTarget = defopts(settings, 'Target', 10^(-8));

  % get function and dimension IDs
  dimInvIds = inverseIndex(funcSet.dims);
  dimIds = dimInvIds(dims);
  funcInvIds = inverseIndex(funcSet.BBfunc);
  funcIds = funcInvIds(BBfunc);

  if ~all(dimIds)
    fprintf('Wrong dimesion request!\n')
  end
  if ~all(funcIds)
    fprintf('Wrong function request!\n')
  end
  
  % count EFE quantiles
  nFunc = length(funcIds);
  nDims = length(dimIds);
  useMaxInstances = 15;
  for alg = 1:numOfData
    for f = 1:nFunc
      for d = 1:nDims
        data_stats{alg}{f, d} = quantile(instanceEFE(data{alg}{funcIds(f), dimIds(d)}, ...
                                         efeTarget, useMaxInstances), quantiles);
      end
    end
  end

  % gain rankings
  nQuantiles = length(quantiles);
  ranks = cell(nFunc, nDims);
  for f = 1:nFunc
    for d = 1:nDims
      % default value is the average rank
      ranks{f,d} = (numOfData+1)/2 * ones(nQuantiles, numOfData);
      notEmptyData = inverseIndex(arrayfun(@(x) ~any(isnan(data_stats{x}{f,d})), 1:numOfData));
      for q = 1:nQuantiles
        thisData = cell2mat(arrayfun(@(x) data_stats{x}{f,d}(q), notEmptyData, 'UniformOutput', false));
        [~, ~, I] = unique(thisData);
        ranks{f,d}(q, notEmptyData) = I';
      end
    end
  end

  % aggregate ranks accross functions
  for dat = 1:numOfData
    for d = 1:nDims
      for ranking = 1:numOfData
        for q = 1:nQuantiles
          rankTable{dat, d}(ranking, q) = sum(arrayfun(@(x) ranks{x,d}(q, dat) == ranking, 1:nFunc));
        end
      end
      rankTable{dat, d} = (1:numOfData)*rankTable{dat,d};
    end
  end
  
  % create sum of ranks table
  rankTable = cell2mat(rankTable);
  sumOfRanks = arrayfun(@(x) sum(rankTable(:, x:nQuantiles:end), 2), 1:nQuantiles, 'UniformOutput', false);
  rankTable = [rankTable, cell2mat(sumOfRanks)];

end

function efe = instanceEFE(fbest, efeTarget, maxInstances)
% Calculates EFE in one instance

  if isempty(fbest)
    efe = NaN;
    return
  end
  
  [maxEvals, nInstances] = size(fbest);
  nInstances = min(nInstances, maxInstances);

  efe = zeros(1, nInstances);
  for i = 1:nInstances
    efe_actual = find(fbest(:, i) < efeTarget, 1, 'first');
    if isempty(efe_actual)
      efe(i) = maxEvals * (1 + 1/9 * log( fbest(end, i) / efeTarget ) );
    else
      efe(i) = efe_actual;
    end
  end
end