classdef PairObliqueSplit < Split
% PairObliqueSplit for each pair of points constructs a normal vector,
% projects all points onto this vector and uses the projected values as
% tresholds for decision to which side the point belongs

  methods
    function obj = PairObliqueSplit(options)
      obj = obj@Split(options);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for i = 1:n-1
        for j = i+1:n
          % normal vector of the hyperplane
          v = obj.X(i, :) - obj.X(j, :);
          % project X onto v
          projectedValues = unique(obj.X * v');
          for treshold = projectedValues
            candidate = obj.splitCandidate;
            candidate.splitter = @(X) ...
              transformApply(X, trans) * v' <= treshold;
            candidate.gain = splitGain.get(candidate.splitter);
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
      end
    end
  end
end