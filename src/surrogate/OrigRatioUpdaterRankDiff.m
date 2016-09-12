classdef OrigRatioUpdaterRankDiff < OrigRatioUpdater
% TODO
% [ ] define different updateRate for positive and negative trend
% [ ] create function errRankMuOnly for two independent f-values vectors
% [ ] weightedSum aggregation of historical errors
    
  properties
    origParams
    lastRatio
    
    parsedParams
    maxRatio
    minRatio
    startRatio
    updateRate
    weights
    ec
    
    rankDiffs
    aggregateType
    lastUpdateGeneration
    
    plotDebug = 0;
    history = [];
    historyRatio = [];
    historyTrend = [];
    fh
  end
  
  methods 
    % get new value of parameter
    function newRatio = update(obj, modelY, origY, ~, ~, countiter, varargin)
      % ratio is updated according to the following formula
      %
      % nHistory = length(weights)
      % newRatio = lastRatio + updateRate*(weights .* rankDiffs((end-nHistory+1):end)     (eqn. 1)
      % newRatio = min(max(newRatio, minRatio), maxRatio)
      %
      % Notes:
      % - update() should be called every generation; if all the samples
      %   were evaluated with the original fitness, model should be also
      %   constructed and reasonable modelY values should be passed here
      % - if update() is not called in any particular generation(s),
      %   it results in NaN entry for that generation(s)
      
      if (nargin >= 7) obj.ec = varargin{1}; end
      obj.rankDiffs((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;
      obj.historyRatio((obj.lastUpdateGeneration+1):(countiter-1)) = obj.lastRatio;
      
      if (isempty(modelY) || (max(modelY) - min(modelY)) == 0 ...
          || (max(origY) - min(origY)) == 0)
        obj.rankDiffs(countiter) = NaN;
      else
        % TODO: create function errRankMuOnly for two independent f-values vectors
        [~, sort1] = sort(modelY);
        ranking2   = ranking(origY);
        rankErr = errRankMuOnly(ranking2(sort1), obj.ec.cmaesState.mu);
        obj.rankDiffs(countiter) = rankErr;
      end
      
      ratio = obj.aggregateTrend();
      
      % obj.lastRatio is initialized as 'startRatio' parameter in the
      % constructor
      % TODO: correct this for faster update!
      lastGenRatio = obj.historyRatio(countiter-1);
      newRatio = (1-obj.updateRate) * lastGenRatio + obj.updateRate * ratio;
      newRatio = min(max(newRatio, obj.minRatio), obj.maxRatio);
      
      if obj.plotDebug
        fprintf('New ratio=%0.2f based on rankDiff trend=%0.2f\n', newRatio, ratio);

        obj.history = [obj.history obj.rankDiffs(countiter)];
        obj.historyTrend = [obj.historyTrend ratio];
      end

      obj.historyRatio(countiter) = newRatio;
      obj.lastRatio = newRatio;
      obj.lastUpdateGeneration = countiter;
    end

    function obj = OrigRatioUpdaterRankDiff(ec, parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      % parameter 'ec' is a reference to the EvolutionControl
      obj.ec = ec;
      obj.parsedParams = struct(parameters{:});
      % maximal possible ratio returned by getValue
      obj.maxRatio = defopts(obj.parsedParams, 'maxRatio', 1);
      % minimal possible ratio returned by getValue
      obj.minRatio = defopts(obj.parsedParams, 'minRatio', 0.02);
      % starting value of ratio for initial generations
      obj.startRatio = defopts(obj.parsedParams, 'startRatio', (obj.maxRatio - obj.minRatio)/2);
      obj.lastRatio  = obj.startRatio;
      % how much is the lastRatio affected by the weighted trend
      obj.updateRate = defopts(obj.parsedParams, 'updateRate', 0.5);
      % type of aggregation of historical values of RankDiff errors
      obj.aggregateType = defopts(obj.parsedParams, 'aggregateType', 'median');
      % weights for weighted sum
      obj.weights = defopts(obj.parsedParams, 'weights', exp([1:4]/2) / sum(exp([1:4]/2)));
      
      obj.rankDiffs = [];
      obj.lastUpdateGeneration = 0;
      
      if obj.plotDebug 
        figure;
        obj.fh = axes;
      end
    end
    
    function value = aggregateTrend(obj)
      % aggregate last criterion values into one value expressing an increasing or
      % decreasing trend
      %
      % This implementation:
      % - takes median of the last length(weights) values
      % - NaN values are ignored (less values are used then)

      nHistory = min(length(obj.rankDiffs), length(obj.weights));
      values = obj.rankDiffs((end-nHistory+1):end);
      bValues  = ~isnan(values);
      if (~any(bValues))
        value = obj.startRatio;
        return;
      end

      switch lower(obj.aggregateType)
      case 'median'
        value = median(values(bValues));
      case 'weightedsum'
        localWeights = obj.weights((end-nHistory+1):end);
        value = sum(localWeights(bValues) .* values(bValues));
      otherwise
        error(sprintf('OrigRatioUpdaterRankDiff: aggregateType ''%s'' not implemented.', obj.aggregateType));
      end
    end
    
    function value = getLastRatio(obj, countiter)
      % TODO: why there is this "+ 1"?!
      if countiter > obj.lastUpdateGeneration + 1
        obj.update([], [], [], [], countiter);
      end
      value = obj.lastRatio;
      
      if obj.plotDebug
          scatter(1:length(obj.history), obj.history, 140, '.');
          hold on;
          scatter(1:length(obj.historyTrend), obj.historyTrend, 140, '.');
          scatter(1:length(obj.historyRatio), obj.historyRatio, 140, '.');
          legend('rankDiff', 'Trend', 'Ratio');
          hold off;
          pause(0.0001);    
      end
    end
    
  end
end
