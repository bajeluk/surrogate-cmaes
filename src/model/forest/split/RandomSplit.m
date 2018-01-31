classdef RandomSplit < Split
% RandomSplit creates a split function used in decision trees randomly
  
  properties %(Access = protected)
    split_nRepeats % number of random repeats
  end
  
  methods
    function obj = RandomSplit(options)
      obj = obj@Split(options);
      obj.split_nRepeats = defopts(options, 'split_nRepeats', 1);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [n, d] = size(obj.split_X);
      for iRepeats = 1:obj.split_nRepeats
        feature = randi(d);
        featureSelector = (1:d == feature)';
        treshold = obj.split_X(randi(n), feature);
        candidate = obj.splitCandidate;
        candidate.splitter = obj.createSplitter(@(X) ...
          X * featureSelector - treshold);
        [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
  
end