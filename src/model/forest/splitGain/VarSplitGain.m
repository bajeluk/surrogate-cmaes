classdef VarSplitGain < SplitGain
% VarSplitGain evaluates split functions used in decision trees using 
% variance of predicted values

  methods
    function obj = VarSplitGain(options)
      obj = obj@SplitGain(options);
      obj.splitGain_computeSd2 = false;
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      value = var(data.y);
    end
  end
end