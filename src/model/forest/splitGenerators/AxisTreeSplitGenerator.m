classdef AxisTreeSplitGenerator < TreeSplitGenerator
  % obj = AxisTreeSplitGenerator(featuresFraction, valuesFraction)
  % for feature = random unique from features
  %   for value = random unique from values(feature)
  %     yield lambda X(:, feature) <= value
  % setting featuresFraction = 1, valuesFraction = 1 generates every
  % possible split point
  properties (Access = protected)
    featuresFraction % fraction of features to split by
    features % sampled features
    iFeatures % current index in features
    valuesFraction % fraction of data values to split by
    values % sampled unique values of current feature
    iValues % current index in values
    pcaTransform
    pcaWeights
  end
  
  methods
    function obj = AxisTreeSplitGenerator(featuresFraction, valuesFraction, pcaTransform)
      if nargin > 0
        obj.featuresFraction = featuresFraction;
        obj.valuesFraction = valuesFraction;
        obj.pcaTransform = pcaTransform;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      if obj.pcaTransform
        obj.pcaWeights = pca(X);
        obj.X = X * obj.pcaWeights;
      else
        obj.X = X;
      end
      obj.y = y;
      obj.features = 1:size(X, 2);
      obj.features = datasample(obj.features, ...
        ceil(obj.featuresFraction * length(obj.features)), ...
        'Replace', false);
      obj.iFeatures = 1;
      obj.iValues = 1;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = obj.iFeatures <= length(obj.features);
    end
    
    function f = next(obj)
      % returns next split function
      feature = obj.features(obj.iFeatures);
      if obj.iValues == 1
        obj.values = unique(obj.X(:, feature));
        obj.values = datasample(obj.values, ...
          ceil(obj.valuesFraction * length(obj.values)), ...
          'Replace', false);
      end
      treshold = obj.values(obj.iValues);
      if obj.pcaTransform
        W = obj.pcaWeights;
        f = @(X) AxisTreeSplitGenerator.selectFeature(X * W, feature) <= treshold;
      else
        f = @(X) AxisTreeSplitGenerator.selectFeature(X, feature) <= treshold;
      end
      obj.iValues = obj.iValues + 1;
      if obj.iValues > length(obj.values)
        obj.iFeatures = obj.iFeatures + 1;
        obj.iValues = 1;
      end
    end
  end
  
  methods (Static, Access = private)
    function x = selectFeature(X, feature)
      x = X(:, feature);
    end
  end
end