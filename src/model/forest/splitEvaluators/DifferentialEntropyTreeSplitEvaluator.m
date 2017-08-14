classdef DifferentialEntropyTreeSplitEvaluator < TreeSplitEvaluator
  methods
    function obj = DifferentialEntropyTreeSplitEvaluator(varargin)
      obj = obj@TreeSplitEvaluator(varargin{:});
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
      value = 0.5 * ( ...
        1 - log(2*pi) + sum(log(data.sd2)) ...
        );
    end
  end
end