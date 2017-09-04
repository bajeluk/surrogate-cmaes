classdef MyHandle < handle
  
  properties (Constant, Access = private)
    nodeTemplate = struct(... % template for nodes
        'X', []);
  end
  
  properties
    X
    sumX
    nodes
    nNodes
  end
  
  methods
    function obj = MyHandle(X)
      obj.X = X;
    end
    
    function s = sumRecursive(obj)
      obj.sumX = 0;
      %obj.nNodes = 1;
      %obj.nodes = repmat(TreeModel.nodeTemplate, 1, 1);
      %obj.nodes(obj.nNodes).X = obj.X;
      
      obj.sumRecursive1(obj.X);
      s = obj.sumX;
    end
    
    function sumRecursive1(obj, X)
      n = size(X, 1);
      if n == 1
        obj.sumX = obj.sumX + X(1, :);
      else
        n1 = floor(n / 2);
        obj.sumRecursive1(X(1:n1, :));
        obj.sumRecursive1(X(n1+1:n, :));
      end
    end
    
    function s = sumIterative(obj)
      obj.sumX = 0;
      pos = zeros(size(obj.X, 1), 1);
      v = (1:size(obj.X, 1))';
      i = 0;
      k = 1;
      while i < k
        idx = pos == i;
        X = obj.X(idx, :);
        n = size(X, 1);
        if n == 1
          obj.sumX = obj.sumX + X(1, :);
        else
          n1 = floor(n / 2);
          idx = v(idx);
          pos(idx(1:n1)) = k;
          pos(idx(n1+1:n)) = k+1;
          k = k+2;
        end
        i = i + 1;
      end
      s = obj.sumX;
    end
  end
  
end

