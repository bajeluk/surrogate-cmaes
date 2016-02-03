classdef OrigRatioUpdaterRMSE < OrigRatioUpdater
  properties
    origParams
    
    parsedParams
    lastRatio
    rmse
    
  end
  
  methods 
    % get new value of parameter
    function newRatio = getValue(obj, modelY, origY, dim, lambda, countiter)
      
      obj.rmse(countiter) = sqrt(sum((modelY' - origY).^2))/length(origY);
      maxR = obj.parsedParams.maxRatio;
      minR = obj.parsedParams.minRatio;
      
      if countiter == 1 || obj.rmse(countiter-1) == 0
        newRatio = obj.parsedParams.startRatio;
      else
        ri = log10(obj.rmse(countiter-1) / obj.rmse(countiter));
        if ri <= -1
          newRatio = maxR;
        elseif ri >= 1
          newRatio = minR;
        else
          newRatio = ri*(minR - maxR)/2 + (minR + maxR)/2;
%         elseif ri > 0
%           newRatio = (1 - ri)*obj.lastRatio;
%         else
%           newRatio = obj.lastRatio*(1 + ri) - ri;
        end
      end
      
      lr = obj.parsedParams.learingRate;
      newRatio = (1-lr)*obj.lastRatio + lr*newRatio;
      obj.lastRatio = newRatio;

    end

    function obj = OrigRatioUpdaterRMSE(parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      obj.parsedParams = struct(parameters{:});
      obj.rmse = [];
    end
  end
end