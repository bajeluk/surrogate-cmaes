classdef DESplitGain < SplitGain
% DESplitGain evaluates split functions used in decision 
% trees using differential entropy

  methods
    function obj = DESplitGain(options)
      obj = obj@SplitGain(options);
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      value = 1 + log(2*pi) + sum(log(data.sd2)) / n;
      value = 0.5 * value;
    end
  end
end