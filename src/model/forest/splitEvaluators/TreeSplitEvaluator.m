classdef (Abstract) TreeSplitEvaluator < handle
  properties (Access = protected)
    X % input data
    y % output data
    gainFunc; % gain function f(current, left, right) -> double
    current % struct holding prediction in current node for gainFunc
  end
  
  methods (Abstract)
    r = eval(f)
    reset(obj, X, y)
    % resets the generator with new data
  end
end