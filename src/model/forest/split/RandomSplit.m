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
      % get number of repeats using only one hyperplane per repetition
      nRepeats = obj.getRepeats(1);
      for iRepeats = 1:nRepeats
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
  
  methods (Access=protected)
    function [nRepeats, maxHypRem] = getRepeats(obj, nHypPerRepeat, maxRepeats)
    % calculate number of repeats and remaining hyperplanes in the last
    % repetition according to the number of hyperplanes per repetition
    
      [n, d] = size(obj.split_X);
      % get the number of hyperplanes
      maxHyp = obj.getMaxHyp(n, d);
      if nargin < 3
        maxRepeats = maxHyp;
      end
      % the number of repeats is limited by its maximal value
      nRepeats = min(obj.split_nRepeats, maxRepeats);
      % gain the number of available repeats
      if nRepeats*nHypPerRepeat > maxHyp
        nRepeats = ceil(maxHyp / nHypPerRepeat);
        % remaining hyperplanes in last repeat
        maxHypRem = mod(maxHyp, nHypPerRepeat);
        if maxHypRem == 0
          maxHypRem = nHypPerRepeat;
        end
      else
        maxHypRem = nHypPerRepeat;
      end
      
    end
  end
  
end