classdef VarSplitGain < SplitGain
% VarSplitGain evaluates split functions used in decision trees using sum
% of variances

  methods
    function obj = VarSplitGain(options)
      obj = obj@SplitGain(options);
      obj.computeSd2 = true;
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      value = sum(data.sd2) / (n^2);
    end
  end
end