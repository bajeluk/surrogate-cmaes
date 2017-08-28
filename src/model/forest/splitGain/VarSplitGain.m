classdef VarSplitGain < SplitGain
% VarSplitGain evaluates split functions used in decision trees using sum
% of variances

  methods
    function obj = VarSplitGain(varargin)
      obj = obj@SplitGain(varargin{:});
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(y);
      value = sum(data.sd2) / (n^2);
    end
  end
end