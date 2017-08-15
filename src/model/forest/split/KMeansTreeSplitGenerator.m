classdef KMeansTreeSplitGenerator < TreeSplitGenerator
  
  properties (Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    metric %
    nRepeats
    iRepeats
    Z % scaled input-output matrix
  end
  
  methods
    function obj = KMeansTreeSplitGenerator(discrimType, metric, nRepeats)
      if nargin > 0
        obj.discrimType = discrimType;
        obj.metric = metric;
        obj.nRepeats = nRepeats;
      end
    end
    
    function reset(obj, X, y)
      % resets the generator with new data
      obj.X = X;
      obj.y = y;
      % gaussians are fit on scaled input-output space
      obj.Z = zscore([X y]);
      obj.iRepeats = 1;
    end
    
    function r = hasNext(obj)
      % whether next split function is available
      r = obj.iRepeats <= obj.nRepeats;
    end
    
    function f = next(obj)
      obj.iRepeats = obj.iRepeats + 1;
      c = kmeans(obj.Z, 2, 'Distance', obj.metric);
      % discriminant analysis of two clusters
      model = fitcdiscr(X, c, 'DiscrimType', obj.discrimType);
      f = @(X) model.predict(X) == 1;
    end
  end
end