classdef RandomSplit < Split
% RandomSplit creates a split function used in decision trees randomly
  
  properties %(Access = protected)
    nRepeats % number of random repeats
  end
  
  methods
    function obj = RandomSplit(options)
      obj = obj@Split(options);
      obj.nRepeats = defopts(options, 'nRepeats', 1);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      [n, d] = size(obj.X);
      for iRepeats = 1:obj.nRepeats
        feature = randi(d);
        featureSelector = (1:d == feature)';
        treshold = obj.X(randi(n), feature);
        candidate = obj.splitCandidate;
        candidate.splitter = obj.createSplitter(@(X) ...
          X * featureSelector - treshold);
        candidate.gain = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
  
end