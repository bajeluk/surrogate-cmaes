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
  end
  
  methods
    function obj = Split(transformationOptions)
      obj.transformationOptions = transformationOptions;
    end
    
    function obj = reset(obj, X, y)
    % sets new transformed input
      [obj.X, obj.y, obj.transformation] = ...
        transform(X, y, obj.transformationOptions);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      candidate = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      feature = randi(d);
      featureSelector = (1:d == feature)';
      value = obj.X(randi(n), feature);
      candidate.splitter = @(X)...
        transformApply(X, trans) * featureSelector <= value;
      candidate.gain = splitGain.get(candidate.splitter);
      best = candidate;
    end
  end
end