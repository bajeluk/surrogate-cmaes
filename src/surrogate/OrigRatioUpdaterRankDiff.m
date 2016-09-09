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
    updateRate
    weights
    ec
    
    rankDiffs
    lastUpdateGeneration
    
    plotDebug = 0;
    historyKendall = [];
    historyRatio = [];
    historyTrend = [];
    fh
  end
  
  methods 
    % get new value of parameter
    function newRatio = update(obj, modelY, origY, ~, ~, countiter)
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
      %   it results in zero NaN entry for that generation(s)
      
      obj.rankDiffs((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;
      
      obj.lastUpdateGeneration = countiter;
            
      if (isempty(modelY) || std(modelY) == 0 || std(origY) == 0)
        obj.rankDiffs(countiter) = NaN;
      else
        % TODO: create function errRankMuOnly for two independent f-values vectors
        [~, sort1] = sort(modelY);
        ranking2   = ranking(origY);
        rankErr = errRankMuOnly(ranking2(sort1), obj.ec.cmaesState.mu);
        obj.rankDiffs(countiter) = rankErr;
      end
      
      ratio = aggregateTrend(obj);
      
      % obj.lastRatio is initialized as 'startRatio' parameter in the
      % constructor
      % TODO: correct this for faster update!
      newRatio = obj.lastRatio + obj.updateRate * ratio;
      newRatio = min(max(newRatio, obj.minRatio), obj.maxRatio);
      
      if obj.plotDebug
          fprintf('New ratio=%0.2f based on kendall trend=%0.2f\n', newRatio, ratio);

          obj.historyRatio = [obj.historyRatio newRatio];
          obj.historyKendall = [obj.historyKendall obj.kendall(countiter)];
          obj.historyTrend = [obj.historyTrend ratio];
      end
      
      obj.lastRatio = newRatio;
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
      obj.lastRatio = defopts(obj.parsedParams, 'startRatio', (obj.maxRatio - obj.minRatio)/2);
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

      switch lower(obj.aggregateType)
      case 'median'
        nHistory = min(length(obj.rankDiffs), length(obj.weights));
        value = median(obj.rankDiffs((end-nHistory+1):end));
      % case 'weightedsum'
      %   % TODO: write this!!
      %   localRankDiffs = ones(1, length(obj.weights));
      %   value = 
      otherwise
        error(sprintf('OrigRatioUpdaterRankDiff: aggregateType ''%s'' not implemented.', obj.aggregateType));
      end
    end
    
    function value = getLastRatio(obj, countiter)
      if countiter > obj.lastUpdateGeneration + 1
        obj.update([], [], [], [], countiter);
      end
      value = obj.lastRatio;
      
      if obj.plotDebug
          scatter(1:length(obj.historyKendall), obj.historyKendall, 140, '.');
          hold on;
          scatter(1:length(obj.historyTrend), obj.historyTrend, 140, '.');
          scatter(1:length(obj.historyRatio), obj.historyRatio, 140, '.');
          legend('Kendall', 'Trend', 'Ratio');
          hold off;
          pause(0.0001);    
      end
    end
    
  end
end
