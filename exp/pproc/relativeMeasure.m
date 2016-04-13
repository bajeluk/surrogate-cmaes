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
  nInstances = defopts(settings, 'MaxInstances', 15);
  refMeasure = defopts(settings, 'RefMeasure', @(x, y) min(x, [], y));
  aggregateDims = defopts(settings, 'AggregateDims', false);
  defaultEvals = defopts(settings, 'DefaultTargets', [2*(1:25), 5*(11:20), 10*(11:25)]);
  
  % count each data median
  data_stats = cellfun(@(D) gainStatistic(D, dimId, funcId, nInstances, false, @median), ...
                            data, 'UniformOutput', false);
  
  % compute reference data
  [nFunc, nDims] = size(data_stats{1});
  for f = 1 : nFunc
    for d = 1 : nDims
      % gather all data from function f and dimension d
      nonEmptyData = ~cellfun(@(D) isempty(D{f, d}), data_stats);
      allActualData = cell2mat(cellfun(@(D) D{f, d}, data_stats(nonEmptyData), 'UniformOutput', false));
      % count reference data
      refData = refMeasure(allActualData, 2);
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
      % bajeluk logarithmic version
%       yTargets = logspace(log(max(refData)), log(max(min(refData), targetValue)), length(defaultEvals));
      % gain target y-values
      if isempty(maxEval)
        yTargets{d} = refData(transformedEvals);  
      else
        yTargets{d} = [refData(transformedEvals(1:end-1)); targetValue];
      end
      
      for D = 1 : length(data_stats)
        targetEvals{D}{f,d} = [];
        if nonEmptyData(D)
          for t = 1 : length(yTargets{d})
            targetEvaluation = find(data_stats{D}{f,d} <= yTargets{d}(t), 1, 'first');
            if isempty(targetEvaluation)
              targetEvals{D}{f,d}(end+1) = NaN;
            else
              targetEvals{D}{f,d}(end+1) = targetEvaluation;
            end
          end
        end
      end
    end
    
    % aggregate results across dimensions
    if aggregateDims
      targetEvals{D}{f} = ceil(funReplaceNaN(@median, cell2mat(targetEvals{D}(f, :)')));
    end
    
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