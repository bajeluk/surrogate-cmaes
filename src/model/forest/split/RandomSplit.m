classdef RandomSplit < Split
% RandomSplit creates a split function used in decision trees randomly
  
  properties %(Access = protected)
    nRepeats % number of random repeats
  end
  
  methods
    function obj = RandomSplit(transformationOptions, nRepeats)
      obj = obj@Split(transformationOptions);
      obj.nRepeats = nRepeats;
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      candidate = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for iRepeats = 1:obj.nRepeats
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
end