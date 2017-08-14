classdef ObliquePairTreeSplitGenerator < TreeSplitGenerator
  % obj = AxisTreeSplitGenerator(featuresFraction, valuesFraction)
  % for feature = random unique from features
  %   for value = random unique from values(feature)
  %     yield lambda X(:, feature) <= value
  % setting featuresFraction = 1, valuesFraction = 1 generates every
  % possible split point
  properties (Access = protected)
    samplesFraction % fraction of data values to split by
    samples % sampled examples
    generated
    evaluator
  end
  
  methods
    function obj = ObliquePairTreeSplitGenerator(evaluator, samplesFraction)
      if nargin > 0
        obj.samplesFraction = samplesFraction;
        obj.evaluator = evaluator;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      obj.X = X;
      obj.y = y;
      obj.samples = datasample(X, ...
        ceil(obj.samplesFraction * size(X, 1)), ...
        'Replace', false);
      obj.generated = false;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = ~obj.generated;
    end
    
    function f = next(obj)
      % returns next split function
      for iSample = 1:size(obj.samples, 1)-1
        for jSample = iSample+1:size(obj.samples, 1)
          u = obj.samples(iSample);
          v = obj.samples(jSample);
          n = u - v;
          projectedValues = unique(obj.samples * n');
          for treshold = projectedValues
            f = @(X) X * n' <= treshold;
            obj.evaluator.eval(f);
          end
        end
      end
      f = obj.evaluator.best.splitter;
    end
  end
end