classdef OrigRatioUpdaterRMSE < OrigRatioUpdater
  properties
    origParams
    lastRatio
    
    parsedParams
    maxRatio
    minRatio
    updateRate
    logRMSEWeights
    logRMSERatioTreshold
    
    rmse
    lastUpdateGeneration
  end
  
  methods 
    % get new value of parameter
    function newRatio = update(obj, modelY, origY, ~, ~, countiter)
      % ratio is updated according to the following formula
      %
      % newRatio = lastRatio + updateRate*(logRMSEWeights.*logRMSE - logRMSERatioTreshold)     (eqn. 1)
      % newRatio = min(max(newRatio, minRatio), maxRatio)
      %
      % Notes:
      % - update() should be called every generation; if all the samples
      %   were evaluated with the original fitness, model should be also
      %   constructed and reasonable modelY values should be passed here
      % - if update() is not called in any particular generation(s),
      %   it results in zero NaN entry for that generation(s)
      
      obj.rmse((obj.lastUpdateGeneration+1):(countiter-1)) = NaN;
      
      obj.lastUpdateGeneration = countiter;
            
      if isempty(modelY)
        obj.rmse(countiter) = NaN;
      else      
        obj.rmse(countiter) = sqrt(sum((modelY - origY).^2))/length(origY);
      end
      
      trend = aggregateRMSETrend(obj);
      
      % obj.lastRatio is initialized as 'startRatio' parameter in the
      % constructor
      newRatio = obj.lastRatio + obj.updateRate * (trend - obj.logRMSERatioTreshold);
      newRatio = min(max(newRatio, obj.minRatio), obj.maxRatio);
      
      obj.lastRatio = newRatio;
    end

    function obj = OrigRatioUpdaterRMSE(parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      obj.parsedParams = struct(parameters{:});
      % maximal possible ratio returned by getValue
      obj.maxRatio = defopts(obj.parsedParams, 'maxRatio', 1);
      % minimal possible ratio returned by getValue
      obj.minRatio = defopts(obj.parsedParams, 'minRatio', 0.1);
      % starting value of ratio for initial generations
      obj.lastRatio = defopts(obj.parsedParams, 'startRatio', (obj.maxRatio - obj.minRatio)/2);
      % how much is the lastRatio affected by the weighted RMSE trend
      obj.updateRate = defopts(obj.parsedParams, 'updateRate', (obj.maxRatio - obj.minRatio)/2);
      % weights for the weighted sum of the log RMSE ratios
      obj.logRMSEWeights = defopts(obj.parsedParams, 'logRMSEWeights', [0.5, 0.3, 0.2]);
      % normalize weights
      obj.logRMSEWeights = obj.logRMSEWeights / sum(obj.logRMSEWeights);
      % minimal value of log RMSE ratio which starts to increase newRatio
      % (update term in (eqn. 1) is positive)
      obj.logRMSERatioTreshold = defopts(obj.parsedParams, 'logRMSERatioTreshold', -0.1);
      obj.rmse = [];
      obj.lastUpdateGeneration = 0;
    end
    
    function value = aggregateRMSETrend(obj)
      % aggregate last RMSE's into one value expressing an increasing or
      % decreasing trend
      %
      % This implementation:
      % - replaces NaN's with maximal values from the last RMSE entries
      % - calculates aggregate value as
      %
      %   partsum(i) = weights(i) * log(rmse(end-i+1) / rmse(end-i)
      %   value      = sum( partsum )
      
      % default value of imaginary RMSE when all RMSE values are NaN
      % in the last obj.rmse entries
      MAX_RMSE = 10^5;
      
      nWeights = min(length(obj.rmse) - 1, length(obj.logRMSEWeights));
      localRMSE = obj.rmse(end - nWeights : end);
      
      local_max_rmse = 2*max(localRMSE(~isnan(localRMSE)));
      if isempty(local_max_rmse)
        localRMSE(isnan(localRMSE)) = MAX_RMSE;
      else
        localRMSE(isnan(localRMSE)) = local_max_rmse;
      end
      
      if length(localRMSE) <= 1
        % we don't have enough RMSE history values ==> stay at the
        % current origRatio ==> set aggregateRMSETrend = 0
        value = 0;
        return
      end
      assert(length(localRMSE)-1 == nWeights, 'DEBUG assertion failed: length of RMSE ~= nWeights + 1');
      
      % replace zeros for division
      localRMSE(localRMSE < eps) = 100*eps;
      
      ratios = log(localRMSE(2:end) ./ localRMSE(1:end-1));
      
      value = sum(obj.logRMSEWeights(1:nWeights) .* ratios(end:-1:1));
    end
    
    function value = getLastRatio(obj, countiter)
      if countiter > obj.lastUpdateGeneration + 1
        obj.update([], [], [], [], countiter)
      end
      value = obj.lastRatio;
    end
    
  end
end
