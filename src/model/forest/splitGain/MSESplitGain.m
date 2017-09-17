classdef MSESplitGain < SplitGain
% MSESplitGain evaluates split functions used in decision trees using MSE

  methods
    function obj = MSESplitGain(options)
      obj = obj@SplitGain(options);
      obj.computeMse = true;
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      %r = data.y - data.yPred;
      %value = r' * r / numel(r);
      value = data.mse;
    end
  end
end