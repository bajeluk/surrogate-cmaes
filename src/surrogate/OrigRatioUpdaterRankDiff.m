classdef OrigRatioUpdaterRankDiff < OrigRatioUpdater
% TODO
% [ ] define different updateRate for positive and negative trend
% [ ] create function errRankMuOnly for two independent f-values vectors
% [ ] weightedSum aggregation of historical errors
% [ ] faster increased updates then decreased (is it really ok?)
    
  properties
    origParams
    lastRatio
    
    surrogateOpts
    ec

    maxRatio
    minRatio
    startRatio
    updateRate
    lowRank
    highRank
    aggregateType
    weights

    rankDiffs
    lastUpdateGeneration
    gain
    newRatio

    plotDebug = 0;
    historyRatios = [];
    historyAggRankDiffs = [];
    fh
  end
  
  methods 
    % get new value of parameter
    function ratio = update(obj, modelY, origY, ~, ~, countiter, varargin)
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
      obj.historyRatios((obj.lastUpdateGeneration+1):(countiter-1)) = obj.lastRatio;
      
      rankErr = NaN;
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
      
      % Decide the best new ratio based on aggregated rankDiff error
      aggRankDiff = obj.aggregateWithHistory();
      if (~isnan(aggRankDiff))
        obj.gain = min(max(0, aggRankDiff - obj.lowRank), (obj.highRank - obj.lowRank)) / (obj.highRank - obj.lowRank);
        obj.newRatio = obj.minRatio + obj.gain * (obj.maxRatio - obj.minRatio);
      else
        obj.gain = NaN;
        obj.newRatio = obj.startRatio;
      end
      % Debug:
      % fprintf('rankErr = %.2f ;  value = %.2f ;  ratio = %.2f\n', rankErr, obj.gain, obj.newRatio);

      % lastGenRatio -- the last used ratio (from the last updated generation)
      if (countiter > 1)
        lastGenRatio = obj.historyRatios(countiter-1);
      else
        lastGenRatio = obj.startRatio;
      end
      % final ratio is exponentially updated:   r = (1-a) * r_old  +  a * r_new
      ratio = (1-obj.updateRate) * lastGenRatio + obj.updateRate * obj.newRatio;
      ratio = min(max(ratio, obj.minRatio), obj.maxRatio);
      
      obj.historyRatios(countiter) = ratio;
      obj.lastRatio = ratio;
      obj.lastUpdateGeneration = countiter;

      if obj.plotDebug
        fprintf('New ratio=%0.2f based on rankDiff trend=%0.2f\n', ratio, obj.newRatio);
        obj.historyAggRankDiffs(countiter) = aggRankDiff;
      end
    end

    function obj = OrigRatioUpdaterRankDiff(ec, parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      % parameter 'ec' is a reference to the EvolutionControl
      obj.ec = ec;
      obj.surrogateOpts = parameters;
      % maximal possible ratio returned by getValue
      obj.maxRatio = defopts(obj.surrogateOpts, 'DTAdaptive_maxRatio', 0.9);
      % minimal possible ratio returned by getValue
      obj.minRatio = defopts(obj.surrogateOpts, 'DTAdaptive_minRatio', 0.02);
      % starting value of ratio for initial generations
      obj.startRatio = defopts(obj.surrogateOpts, 'DTAdaptive_startRatio', (obj.maxRatio - obj.minRatio)/2);
      obj.lastRatio  = obj.startRatio;
      % how much is the lastRatio affected by the weighted trend
      obj.updateRate = defopts(obj.surrogateOpts, 'DTAdaptive_updateRate', 0.9);
      % type of aggregation of historical values of RankDiff errors
      obj.aggregateType = defopts(obj.surrogateOpts, 'DTAdaptive_aggregateType', 'weightedSum');
      % weights for weighted sum
      obj.weights = defopts(obj.surrogateOpts, 'DTAdaptive_weights', exp([1:4]/2) / sum(exp([1:4]/2)));
      % lowest and highest rank which affect gain util it saturates to 0 or 1
      obj.lowRank  = defopts(obj.surrogateOpts, 'DTAdaptive_lowRank', 0.1);
      obj.highRank = defopts(obj.surrogateOpts, 'DTAdaptive_highRank', 0.5);
      
      obj.rankDiffs = [];
      obj.lastUpdateGeneration = 0;
      
      if obj.plotDebug 
        figure;
        obj.fh = axes;
      end
    end
    
    function value = aggregateWithHistory(obj)
      % aggregate last criterion values into one value
      %
      % This implementation:
      % - ignores NaN values in history (less values are used then)
      %
      % (a) takes median of the last length(weights) values
      %     or
      % (b) takes weighted sum of the last length(weights) values

      % take at most nHistory last values
      nHistory = min(length(obj.rankDiffs), length(obj.weights));
      values = obj.rankDiffs((end-nHistory+1):end);
      % identify NaN's
      bValues  = ~isnan(values);
      % return with NaN if no valid values in history
      if (~any(bValues))
        value = NaN;
        return;
      end

      switch lower(obj.aggregateType)
      case 'median'
        value = median(values(bValues));
      case 'weightedsum'
        % take adequate weights and re-norm them to sum to 1
        localWeights = obj.weights((end-nHistory+1):end);
        localWeights = localWeights(bValues) ./ sum(localWeights(bValues));
        % return the weighted sum
        value = sum(localWeights .* values(bValues));
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
          scatter(1:length(obj.historyAggRankDiffs), obj.historyAggRankDiffs, 140, '.');
          scatter(1:length(obj.historyRatios), obj.historyRatios, 140, '.');
          legend('rankDiff', 'Trend', 'Ratio');
          hold off;
          pause(0.0001);    
      end
    end
    
  end
end
