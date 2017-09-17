classdef DEMSDSplitGain < SplitGain
% UMVUESplitGain evaluates split functions used in decision trees using
% uniformly minimum-variance unbiased estimator

  methods
    function obj = DEMSDSplitGain(options)
      obj = obj@SplitGain(options);
      obj.computeSd2 = true;
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      value = log(exp(1)*pi) ...
        + log(sum(data.sd2) + realmin) ...
        - digamma((n - 1) / 2);
      value = 0.5 * value;
    end
  end
end