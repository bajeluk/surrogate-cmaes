classdef Value
  properties
    X
  end
  
  methods
    function obj = Value(X)
      obj.X = X;
    end
    
    function obj = set(obj, i, j, k)
      obj.X(i, j) = k;
    end
  end
  
end