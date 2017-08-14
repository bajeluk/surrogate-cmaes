classdef AGEntropyTreeSplitEvaluator < TreeSplitEvaluator
  methods
    function obj = AGEntropyTreeSplitEvaluator()
      obj = obj@TreeSplitEvaluator();
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
      value = 0.5 * (...
        log(exp(1)*pi) + log(sum(data.y.^2)) + digamma(n/2)...
        ); 
    end
  end
end