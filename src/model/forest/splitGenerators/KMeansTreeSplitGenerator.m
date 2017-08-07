classdef KMeansTreeSplitGenerator < TreeSplitGenerator
  % obj = KMeansTreeSplitGenerator(featuresFraction, valuesFraction)
  % for feature = random unique from features
  %   for value = random unique from values(feature)
  %     yield lambda X(:, feature) <= value
  % setting featuresFraction = 1, valuesFraction = 1 generates every
  % possible split point
  properties (Access = protected)
    metric
    normalize
    xMean
    xVar
    nRepeats
    iRepeats
  end
  
  methods
    function obj = KMeansTreeSplitGenerator(metric, normalize, nRepeats)
      if nargin > 0
        obj.metric = metric;
        obj.normalize = normalize;
        obj.nRepeats = nRepeats;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      if obj.normalize
        obj.xMean = mean(X);
        obj.xVar = var(X) + 1e-10;
        obj.X = bsxfun(@rdivide, bsxfun(@minus, X, obj.xMean), obj.xVar);
      else
        obj.X = X;
      end
      obj.y = y;
      obj.iRepeats = 1;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = obj.iRepeats <= obj.nRepeats;
    end
    
    function f = next(obj)
      obj.iRepeats = obj.iRepeats + 1;
      [~, C] = kmeans(obj.X, 2, 'Distance', obj.metric);
      metric = obj.metric;
      if strcmpi(metric, 'sqeuclidean')
        metric = 'squaredeuclidean';
      end
      if obj.normalize
        xMean = obj.xMean;
        xVar = obj.xVar;
        f = @(X) KMeansTreeSplitGenerator.split(bsxfun(@rdivide, bsxfun(@minus, X, xMean), xVar), C, metric);
      else
        f = @(X) KMeansTreeSplitGenerator.split(X, C, metric);
      end
    end
  end
  
  methods (Static, Access = private)
    function P = split(X, C, metric)
      % returns 1 if X(i,:) is closer to C(1,:) than to C(2,:)
      D = pdist2(X, C, metric);
      P = D(:, 1) < D(:, 2);
    end
  end
end