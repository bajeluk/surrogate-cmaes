classdef DESplitGain < SplitGain
% DESplitGain evaluates split functions used in decision 
% trees using differential entropy

  methods
    function obj = DESplitGain(options)
      obj = obj@SplitGain(options);
      obj.splitGain_computeSd2 = true;
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      value = 1 + log(2*pi) + sum(log(data.sd2 + realmin)) / (2*n);
    end
  end
end