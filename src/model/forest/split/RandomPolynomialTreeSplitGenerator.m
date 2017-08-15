classdef RandomPolynomialTreeSplitGenerator < TreeSplitGenerator
  
  properties
    nFeaturesInteract
    nRepeats
    iRepeats
    degree
  end
  
  methods
    function obj = RandomPolynomialTreeSplitGenerator(...
        degree, nFeaturesInteract, nRepeats)
      if nargin > 0
        obj.nFeaturesInteract = nFeaturesInteract;
        obj.nRepeats = nRepeats;
        obj.degree = degree;
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
      switch obj.degree
        case 'linear2'
          % each column has values within [min(column), max(column)]
          featuresRand = rand(nFeatures, nFeatures) ...
            .* repmat(featuresMax - featuresMin, nFeatures, 1) ...
            + repmat(featuresMin, nFeatures, 1);
          line = featuresRand \ ones(nFeatures, 1);
          f = @(X) X(:, features) * line <= 1;
        otherwise
          % select point where hyperplane passes through
          origin = rand(1, nFeatures) ...
            .* (featuresMax - featuresMin) ...
            + featuresMin;
          %origin = datasample(X, 1);
          degree = obj.degree;
          nFeaturesP = size(generateFeatures(features, degree, false), 2);
          % select a direction of the hyperplane
          angles = rand(nFeaturesP, 1) * pi - pi/2;
          % convert direction to weights
          weights = tan(angles);
          f = @(X) ...
            generateFeatures(bsxfun(@minus, X(:, features), origin), degree, false)...
            * weights <= 0;
      end
    end
  end
end

