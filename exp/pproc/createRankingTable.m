function [rankTable, ranks] = createRankingTable(data, varargin)
% [table, ranks] = createRankingTable(data, settings)
% Creates table containing rankings for different evaluations.
%
% Input:
%   data      - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%     'Evaluations' - evaluations chosen to count
%     'Statistic'   - statistic of data | string or handle (@mean, 
%                       @median)
%
% Output:
%   table - table of rankings
%   ranks - rankings for each function and dimension
%
% See Also:
%   rankingTable, speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  rankTable = [];
  if nargin < 1 || isempty(data)
    help createRankingTable
    return
  end
  settings = settings2struct(varargin{:});
  
  numOfData = length(data);
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  dims    = defopts(settings, 'TableDims', funcSet.dims);
  BBfunc  = defopts(settings, 'TableFuns', funcSet.BBfunc);
  evaluations = defopts(settings, 'Evaluations', [20 40 80]);
  statistic = defopts(settings, 'Statistic', @median);
  if ischar(statistic)
    if strcmp(statistic, 'quantile')
      statistic = @(x, dim) quantile(x, [0.25, 0.5, 0.75], dim);
    else
      statistic = str2func(statistic);
    end
  end

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

  % count means
  useMaxInstances = 15;
  data_stats = cellfun(@(D) gainStatistic(D, dimIds, funcIds, ...
                            'MaxInstances', useMaxInstances, ...
                            'AverageDims', false, ...
                            'Statistic', statistic, ...
                            'SuppWarning', true), ...
                       data, 'UniformOutput', false);
                     
  % gain rankings
  nFunc = length(funcIds);
  nDims = length(dimIds);
  nEvals = length(evaluations);
  ranks = cell(nFunc, nDims);
  for f = 1:nFunc
    for d = 1:nDims
      % init
      ranks{f,d} = zeros(nEvals, numOfData);
      notEmptyData = inverseIndex(arrayfun(@(x) ~isempty(data_stats{x}{f,d}), 1:numOfData));
      for e = 1:nEvals
        thisData = cell2mat(arrayfun(@(x) data_stats{x}{f,d}(evaluations(e)), notEmptyData, 'UniformOutput', false));
        thisData = max(thisData, 1e-8 * ones(size(thisData)));
        [~, ~, I] = unique(thisData);
        ranks{f,d}(e, notEmptyData) = I';
      end
    end
  end

  % aggregate ranks accross functions
  for dat = 1:numOfData
    for d = 1:nDims
      for ranking = 1:numOfData
        for e = 1:nEvals
          rankTable{dat, d}(ranking, e) = sum(arrayfun(@(x) ranks{x,d}(e, dat) == ranking, 1:nFunc));
        end
      end
    end
  end
  
  % create first rank table
  rankTable = createTable(rankTable, 1);

end

function rankTable = createTable(table, rank)
% Creates table of sums of rank
  
  [numOfData, nDims] = size(table);
  nEvals = size(table{1,1}, 2);
  
  rankTable = zeros(numOfData, (nDims+1)*nEvals);
  % data
  for dat = 1:numOfData
    % dimensions
    for d = 1:nDims
      % evaluations
      for e = 1:nEvals
        % gain ranks
        rankTable(dat, (d-1)*nEvals + e) = table{dat,d}(rank, e);
      end
    end
    % rank sums
    for e = 1:nEvals
      % print only first ranks
      rankTable(dat, nDims*nEvals + e) = sum(arrayfun(@(x) table{dat,x}(rank, e), 1:nDims));
    end
  end
end