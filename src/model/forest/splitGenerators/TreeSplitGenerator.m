classdef (Abstract) TreeSplitGenerator < Iterator
  properties (Access = protected)
    X % input data
    y % output data
  end
  
  methods (Abstract)
    reset(obj, X, y)
    % resets the generator with new data
  end
end