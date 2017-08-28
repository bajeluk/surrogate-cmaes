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
      n = size(data.y, 1);
      value = 0.5 * (...
        log(exp(1)*pi) + log(sum(data.y.^2)) - digamma(n/2)...
        ); 
    end
  end
end