classdef SSESplitGain < SplitGain
% SSESplitGain evaluates split functions used in decision trees using SSE

  methods
    function obj = SSESplitGain(options)
      obj = obj@SplitGain(options);
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      r = data.y - data.yPred;
      value = r' * r / numel(r);
    end
  end
end