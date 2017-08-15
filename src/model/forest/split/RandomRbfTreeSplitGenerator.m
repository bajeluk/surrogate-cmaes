classdef RandomRbfTreeSplitGenerator < TreeSplitGenerator
  
  properties
    nFeaturesInteract
    nRepeats
    iRepeats
    metric
  end
  
  methods
    function obj = RandomRbfTreeSplitGenerator(...
        metric, nFeaturesInteract, nRepeats)
      if nargin > 0
        obj.nFeaturesInteract = nFeaturesInteract;
        obj.nRepeats = nRepeats;
        obj.metric = metric;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      obj.X = X;
      obj.y = y;
      obj.iRepeats = 1;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = obj.iRepeats <= obj.nRepeats ...
        && obj.nFeaturesInteract > 0;
    end
    
    function f = next(obj)
      obj.iRepeats = obj.iRepeats + 1;
      % returns next split function
      nFeatures = min(obj.nFeaturesInteract, size(obj.X, 2));
      features = datasample(1:size(obj.X, 2), nFeatures, 'Replace', false);
      X = obj.X(:, features);
      featuresMin = min(X);
      featuresMax = max(X);
      % select point where the RBF oiginates
      origin = rand(1, nFeatures) ...
        .* (featuresMax - featuresMin) ...
        + featuresMin;
      %origin = datasample(X, 1);
      metric = obj.metric;
      distances = pdist2(X, origin, metric);
      maxDistance = max(distances);
      minDistance = min(distances);
      treshold = rand() * (maxDistance - minDistance) + minDistance;
      switch metric
        case 'mahalanobis'
          C = nancov(X);
          f = @(X) pdist2(X(:, features), origin, metric, C) <= treshold;
        otherwise
          f = @(X) pdist2(X(:, features), origin, metric) <= treshold;
      end
    end
  end
end

