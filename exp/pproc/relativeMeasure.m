function [targetEvals, yTargets] = relativeMeasure(data, dimId, funcId, defaultEvals, averageDims, refMeasure, nInstances, targetValue)
% Returns cell array of relative statistics accross chosen dimensions 
% for each function
%
% Input:
%   refMeasure - handle to measure of reference data

  if nargin < 8 
    targetValue = 10^-8;
    if nargin < 7
      nInstances = 15;
      if nargin < 6
        refMeasure = @(x, y) min(x, [], y);
        if nargin < 5
          averageDims = false;
          if nargin < 4
            defaultEvals = [2*(1:25), 5*(11:20), 10*(11:25)];
            if nargin < 1
              help relativeMeasure
              return
            end
          end
        end
      end
    end
  end
  
  % count each data median
  data_stats = cellfun(@(D) gainStatistic(D, dimId, funcId, nInstances, averageDims, @median), ...
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
      % no algorithm reached targetValue
      if isempty(maxEval)
        transformedEvals = defaultEvals;
      % at least one algorithm reached targetValue
      else
        transformedEvals = ceil(defaultEvals/max(defaultEvals) * maxEval);
      end
      % Lukas logarithmic version
%       yTargets = logspace(log(max(refData)), log(max(min(refData), targetValue)), length(defaultEvals));
      % gain target y-values
      yTargets = refData(transformedEvals);
      for D = 1 : length(data_stats)
        targetEvals{D}{f,d} = [];
        if nonEmptyData(D)
          for t = 1 : length(yTargets)
            targetEvaluation = find(data_stats{D}{f,d} <= yTargets(t), 1, 'first');
            if isempty(targetEvaluation)
              targetEvals{D}{f,d}(end+1) = NaN;
            else
              targetEvals{D}{f,d}(end+1) = targetEvaluation;
            end
          end
        end
      end
    end
  end
  
end