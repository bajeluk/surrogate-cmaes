function [rankTable, ranks, values] = createRankingTable(data, varargin)
% [table, ranks, values] = createRankingTable(data, settings)
% Creates table containing rankings for different evaluations.
%
% Input:
%   data     - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'Evaluations' - evaluations chosen to count
%     'Rank'        - rank of resulting table or 'sum' to sum all ranks
%     'Ranking'     - type of ranking (see example below)
%                       'tolerant' - equal rank independence
%                       'precise'  - equal ranks shift following ranks
%                       'median'   - equal ranks replaced by medians of
%                                    shifted ranks (from 'precise')
%     'Statistic'   - statistic of data | string or handle (@mean, 
%                       @median)
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%     'Target'      - target value to reach
%
% Output:
%   table  - table of rankings
%   ranks  - rankings for each function and dimension
%   values - values of data in demanded evaluations
%
% Ranking Example:
%    values  =    [ 1    5    8   13    1    8    1   21 ]
%
%   'tolerant':   [ 1    2    3    4    1    3    1    5 ]
%   'precise':    [ 1    4    5    7    1    5    1    8 ]
%   'median':     [ 2    4   5.5   7    2   5.5   2    8 ]
%
% See Also:
%   rankingTable, createEFETable, speedUpPlot, speedUpPlotCompare, 
%   dataReady

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
  tableMode = defopts(settings, 'Mode', 'evaluations');
  if strcmp(tableMode, 'evaluations')
    evaluations = defopts(settings, 'Evaluations', [20 40 80]);
    if any(evaluations < 1)
      warning(['In table mode ''evaluations'' all evaluation values has to be >= 1.',...
               'Dividing by maximal element (%d)'], max(evaluations))
      evaluations = evaluations/max(evaluations);
    end
  else
    evaluations = defopts(settings, 'Evaluations', [1/3 1]);
  end
  targetValue = defopts(settings, 'Target', 1e-8);
  tableRank = defopts(settings, 'Rank', 1);
  chosenRank = chooseRanking(defopts(settings, 'Ranking', 'tolerant'));
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
  values = cell(nFunc, nDims);
  
  % return ranks in specific evaluations
  for f = 1:nFunc
    for d = 1:nDims
      % init
      ranks{f,d} = zeros(nEvals, numOfData);
      values{f,d} = NaN(nEvals, numOfData);
      % find not empty data
      notEmptyData = inverseIndex(arrayfun(@(x) ~isempty(data_stats{x}{f,d}), 1:numOfData));
      % find maximum of evaluations or return the greatest user-defined 
      % number of evaluations
      maxEvals = max([evaluations(end), arrayfun(@(x) numel(data_stats{x}{f, d}), notEmptyData)]);
      actualData = [];
      % gather actual (function, dimension) non-empty data and fill the
      % missing values to maxEvals with the best achieved result
      for dat = 1:numel(notEmptyData)
        neDat = notEmptyData(dat);
        ds_act = data_stats{neDat}{f, d};
        actualData(1:maxEvals, dat) = [ds_act; ds_act(end)*ones(maxEvals - numel(ds_act), 1)];
      end
%       actualData = cell2mat(arrayfun(@(x) data_stats{x}{f,d}, notEmptyData, 'UniformOutput', false));
      % return ranks in specific evaluations
      if strcmp(tableMode, 'evaluations')
        if all(evaluations <= 1)
          actualEvals = ceil(evaluations*size(actualData, 1));
        else
          actualEvals = evaluations;
        end
      % return ranks when reaching specified target
      else
        % number of evaluations needed by the best data to reach the
        % optimum
        minEval = min([size(actualData, 1), size(actualData, 1) - sum(actualData < targetValue) + 1]);
        actualEvals = ceil(evaluations*minEval);
      end
      
      % asociate ranks to data in specified evaluation numbers
      for e = 1:nEvals
%         thisData = cell2mat(arrayfun(@(x) data_stats{x}{f,d}(actualEvals(e)), notEmptyData, 'UniformOutput', false));
        thisData = actualData(actualEvals(e), :);
        thisData = max(thisData, targetValue * ones(size(thisData)));
        ranks{f,d}(e, notEmptyData) = chosenRank(thisData);
        values{f,d}(e, notEmptyData) = thisData;
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
  
  % create rank table
  rankTable = createTable(rankTable, tableRank, nFunc);

end

function rankTable = createTable(table, rank, nFunc)
% Creates table of sums of rank
  
  [numOfData, nDims] = size(table);
  nEvals = size(table{1,1}, 2);
  sumRank = strcmp(rank, 'sum');
  avgRank = (numOfData + 1)/2;
  
  rankTable = zeros(numOfData, (nDims+1)*nEvals);
  % data
  for dat = 1:numOfData
    % dimensions
    for d = 1:nDims      
      % gain sum of all ranks
      if sumRank
        % find number of missing functions
        missingFuns = nFunc - sum(table{dat,d});
        % count sum of ranks and replace missing functions by the average rank
        rankTable(dat, (d-1)*nEvals+1 : d*nEvals) = (1:numOfData)*table{dat,d} + avgRank*missingFuns;
      % gain chosen ranks
      else
        % evaluations
        for e = 1:nEvals
          rankTable(dat, (d-1)*nEvals + e) = table{dat,d}(rank, e);
        end
      end
    end
    % gain sum of sums of all ranks
    if sumRank
      rankTable(dat, nDims*nEvals + 1: end) = arrayfun(@(x) sum(rankTable(dat, x:nEvals:end)), 1:nEvals);
    % rank sums of chosen ranks
    else
      for e = 1:nEvals
        % print only chosen ranks
        rankTable(dat, nDims*nEvals + e) = sum(arrayfun(@(x) table{dat,x}(rank, e), 1:nDims));
      end
    end
  end
end

function chosenRank = chooseRanking(userRanking)
% Returns appropriate existing ranking to the chosen one

  switch userRanking
    case 'precise'
      chosenRank = @preciseRank; 
    case 'median'
      chosenRank = @medianRank;
    otherwise
      chosenRank = @tolerantRank;
  end

end

function tr = tolerantRank(values)
% Computes ranking of vector elements, where number of equal ranks does not
% play role.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   tr = tolerantRank(a);
%   tr = [ 1    2    3    4    1    3    1    5 ]

  [~, ~, tr] = unique(values);
  tr = tr';

end

function pr = preciseRank(values)
% Computes ranking of vector elements, where number of equal ranks shifts
% the following ranks.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   pr = preciseRank(a);
%   pr = [ 1    4    5    7    1    5    1    8 ]

  tr = tolerantRank(values);
  pr = arrayfun(@(x) sum(tr < x), tr) + 1;

end

function mr = medianRank(values)
% Computes ranking of vector elements, where same ranks are replaced by
% medians of shifted ranks (from preciseRank).
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   mr = medianRank(a);
%   mr = [ 2    4   5.5   7    2   5.5   2    8 ]

  mr = preciseRank(values);
  ranks = unique(mr);
  % equal ranks replace by medians
  if numel(ranks) < numel(mr)
    for r = ranks
      % average rank + actual rank - 1
      mr(mr == r) = (sum(mr == r) + 1)/2 + r - 1;
    end
  end

end