classdef OrigRatioUpdaterRMSE < OrigRatioUpdater
  properties
    origParams
    
    parsedParams
  end
  
  methods 
    % get new value of parameter
    function newRatio = getValue(obj, modelY, origY, dim, lambda, countiter)
      
    end

    function obj = OrigRatioUpdaterRMSE(parameters)
      % constructor
      obj = obj@OrigRatioUpdater(parameters);
      obj.parsedParams = struct(parameters{:});
    end
  end
end