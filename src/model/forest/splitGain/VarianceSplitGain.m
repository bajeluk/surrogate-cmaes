classdef VarianceSplitGain < SplitGain
% SSESplitGain evaluates split functions used in decision trees using sum
% of variances

  methods
    function obj = VarianceSplitGain(varargin)
      obj = obj@SplitGain(varargin{:});
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      value = sum(data.sd2);
    end
  end
end