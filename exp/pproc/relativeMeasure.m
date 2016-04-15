function [targetEvals, yTargets] = relativeMeasure(data, dimId, funcId, varargin)
% Returns cell array of relative statistics accross chosen dimensions 
% for each function
%
% Input:
%   data
%   dimId
%   funId
%   settings:
%     RefMeasure - handle to measure of reference data
%     DefaultTargets

  % initialization
  targetEvals = {};
  yTargets = {};
  if nargin < 1 || isempty(data)
    help relativeMeasure
    return
  end
  if isstruct(varargin)
    settings = varargin;
  else
    % keep cells as cells due to struct command
    vararCellId = cellfun(@iscell, varargin);
    varargin(vararCellId) = {varargin(vararCellId)};
    settings = struct(varargin{:});
  end

  % initialize settings
  targetValue = defopts(settings, 'TargetValue', 10^-8);
  maxInstances = defopts(settings, 'MaxInstances', 15);
  refMeasure = defopts(settings, 'RefMeasure', @(x, y) min(x, [], y));
  aggregateDims = defopts(settings, 'AggregateDims', false);
  defaultEvals = defopts(settings, 'DefaultTargets', [2*(1:25), 5*(11:20), 10*(11:25)]);
  
  nFunc = length(funcId);
  nDims = length(dimId);
  for f = 1 : nFunc
    for d = 1 : nDims
      fId = funcId(f);
      dId = dimId(d);
      % gather all data from function f and dimension d
      nonEmptyData = ~cellfun(@(D) isempty(D{fId, dId}), data);
      allActualData = cellfun(@(D) D{fId, dId}(:, 1:maxInstances), data(nonEmptyData), 'UniformOutput', false);
      % count each data median
      data_stats = cellfun(@(D) funIgnoreNaN(@(x) median(x, 2), D), allActualData, 'UniformOutput', false);
      data_stats = cell2mat(data_stats')';
      % count reference data
      refData = refMeasure(data_stats, 2); 
      % count reference evaluations
      maxEval = find(refData <= targetValue, 1, 'first');
      minEval = find(~isnan(refData), 1, 'first');
      % no algorithm reached targetValue
      if isempty(maxEval) && min(defaultEvals) >= minEval
        transformedEvals = defaultEvals;
      % at least one algorithm reached targetValue
      else
        transformedEvals = minEval + round((defaultEvals - min(defaultEvals))/(max(defaultEvals) - min(defaultEvals)) * (maxEval - minEval));
      end
      
      % gain target y-values
      if isempty(maxEval)
        yTargets{f,d} = refData(transformedEvals);  
      else
        yTargets{f,d} = [refData(transformedEvals(1:end-1)); targetValue];
      end
      
      for D = 1 : length(data)
        [nEvals, nInstances] = size(data{D}{fId, dId});
        nInstances = min(maxInstances, nInstances);
        targetEvals{D}{f,d} = zeros(nInstances, nEvals);
        if nonEmptyData(D)
          for instId = 1:nInstances
            tId = 1;
            dataColumn = data{D}{fId, dId}(:, instId);
            for evalId = 1:nEvals
              if (tId <= length(yTargets{f,d})) && (dataColumn(evalId) <= yTargets{f,d}(tId))
                tId = findNextTarget(yTargets{f,d}, dataColumn(evalId), tId);
              end
              targetEvals{D}{f,d}(instId, evalId) = tId - 1;
            end
          end
        end
        
        % aggregate in one function one dimension
        if ~aggregateDims
          targetEvals{D}{f,d} = ceil(median(targetEvals{D}{f,d}));
        end
      end
    end
    
    % aggregate results across dimensions
    if aggregateDims
      for D = 1 : length(data)
        targetEvals{D}{f} = ceil(median(cell2mat(targetEvals{D}(f, :)')));
      end
    end
    
  end
  
  if aggregateDims
    for D = 1 : length(data)
      targetEvals{D} = targetEvals{D}(:,1);
    end
  end
  
end

function t = findNextTarget(targets, value, tId)
% find next target in 'targets' such that it is <= value starting from 
% index tId
  t = length(targets) + 1;
  for i = tId : length(targets)
    if targets(i) < value
      t = i;
      return
    end
  end
end

function res = funIgnoreNaN(fun, X)
% calculates function ignoring NaN values
% returns NaN only if all values are NaN
  for i = 1:size(X, 1)
    res(i) = fun(X(i, ~isnan(X(i, :))));
  end
end

function res = funReplaceNaN(fun, X)
% calculates function replacing NaN values by maximal values
  maxValue = max(max(X));
  for i = 1:size(X, 2)
    X(isnan(X(:, i)), i) = maxValue;
    res(i) = fun(X(:, i));
  end
end