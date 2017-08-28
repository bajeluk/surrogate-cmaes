classdef UMVUESplitGain < SplitGain
% UMVUESplitGain evaluates split functions used in decision trees using
% uniformly minimum-variance unbiased estimator

  methods
    function obj = UMVUESplitGain()
      obj = obj@SplitGain();
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      value = log(exp(1)*pi) ...
        + log(data.y' * data.y) ...
        - digamma((n + 1 - 1) / 2);
      value = 0.5 * value;
    end
  end
end