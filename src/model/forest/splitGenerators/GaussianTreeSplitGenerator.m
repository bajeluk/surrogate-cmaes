classdef GaussianTreeSplitGenerator < TreeSplitGenerator
  
  properties (Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    nRepeats % number or repeats
    iRepeats % current repeat
    Z % scaled input-output matrix
  end
  
  methods
    function obj = GaussianTreeSplitGenerator(discrimType, nRepeats)
      if nargin > 0
        obj.degree = discrimType;
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
      % fit 2 gaussian distributions on input-output space
      model = fitgmdist(obj.Z, 2, 'RegularizationValue', 0.001);
      c = model.cluster(X);
      % discriminant analysis of two clusters
      model = fitcdiscr(X, c, 'DiscrimType', obj.discrimType);
      f = @(X) model.predict(X) == 1;
    end
  end
end