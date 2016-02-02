classdef OrigRatioUpdater < handle
  properties
    origParams
  end
  
  methods (Abstract)
    % get new value of parameter
    newRatio = getValue(obj, modelY, origY, dim, lambda, countiter);
  end
  
  methods
    function obj = OrigRatioUpdater(parameters)
      % constructor
      obj.origParams = parameters;
    end
  end
end