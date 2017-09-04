classdef Split
% Split creates a split function used in decision trees
  
  properties (Constant)
    splitCandidate = struct( ...
      'splitter', @(X) false(size(X, 1), 1), ...
      'gain', -inf ...
      );
  end
  
  properties %(Access = protected)
    transformationOptions % transform options
    transformation % transformation
    X % input data
    y % output data
    allEqual % whether all y values are equal
  end
  
  methods
    function obj = Split(options)
      obj.transformationOptions = defopts(options, 'transformationOptions', struct);
    end
    
    function obj = reset(obj, X, y)
    % sets new transformed input
      [obj.X, obj.y, obj.transformation] = ...
        transform(X, y, obj.transformationOptions);
      n = size(X, 1);
      varY = var(y);
      obj.allEqual = any(isnan(varY)) || all(varY < eps(max(varY)) * n);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      candidate = obj.splitCandidate;
      candidate.gain = splitGain(candidate.splitter);
      best = candidate;
    end
  end
end