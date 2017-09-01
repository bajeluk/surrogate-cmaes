classdef AxisSplit < Split
% AxisSplit finds the best axis parallel split
  
  methods
    function obj = AxisSplit(options)
      obj = obj@Split(options);
    end

    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for feature = 1:d
        featureSelector = (1:d == feature)';
        values = unique(obj.X(:, feature))';
        for treshold = values
          candidate = obj.splitCandidate;
          candidate.splitter = @(X)...
            transformApply(X, trans) * featureSelector <= treshold;
          candidate.gain = splitGain.get(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
    end
  end
  
end