classdef GMMTreeSplitGenerator < TreeSplitGenerator
  
  properties (Access = protected)
    nRepeats
    iRepeats
    features
  end
  
  methods
    function obj = GMMTreeSplitGenerator(nRepeats)
      if nargin > 0
        obj.nRepeats = nRepeats;
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
      r = obj.iRepeats <= obj.nRepeats;
    end
    
    function f = next(obj)
      obj.iRepeats = obj.iRepeats + 1;
      model = fitgmdist(obj.X, 2, 'RegularizationValue', 0.1);
      f = @(X) model.cluster(X) == 1;
    end
  end
end