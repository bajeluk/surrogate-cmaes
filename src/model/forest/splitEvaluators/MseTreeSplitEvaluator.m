classdef MseTreeSplitEvaluator < TreeSplitEvaluator  
  methods
    function obj = MseTreeSplitEvaluator(varargin)
      obj = obj@TreeSplitEvaluator(varargin{:});
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
      value = sum((data.y - data.yPred).^2);
    end
  end
end