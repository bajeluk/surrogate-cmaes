classdef VarianceTreeSplitEvaluator < TreeSplitEvaluator  
  methods
    function obj = VarianceTreeSplitEvaluator(varargin)
      obj = obj@TreeSplitEvaluator(varargin{:});
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
      value = sum(data.sd2);
    end
  end
end