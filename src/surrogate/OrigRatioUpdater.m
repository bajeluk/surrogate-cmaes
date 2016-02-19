classdef (Abstract) OrigRatioUpdater < handle
  properties (Abstract)
    origParams
    lastRatio
  end
  
  methods (Abstract)
    % get new value of parameter
    newRatio = update(obj, modelY, origY, dim, lambda, countiter);
  end
  
  methods
    function obj = OrigRatioUpdater(parameters)
      % constructor
      obj.origParams = parameters;
    end
    
    function value = getLastRatio(obj, countiter)
      value = obj.lastRatio;
    end
  end
end
