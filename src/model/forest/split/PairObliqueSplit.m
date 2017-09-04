classdef PairObliqueSplit < Split
% PairObliqueSplit for each pair of points constructs a normal vector,
% projects all points onto this vector and uses the projected values as
% tresholds for decision to which side the point belongs

  properties %(Access = protected)
    nQuantize % quantization of tresholds
  end

  methods
    function obj = PairObliqueSplit(options)
      obj = obj@Split(options);
      obj.nQuantize = defopts(options, 'nQuantize', 0);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for i = 1:n-1
        for j = i+1:n
          % normal vector of the hyperplane
          v = obj.X(i, :) - obj.X(j, :);
          % project X onto v
          values = (obj.X * v')';
          if obj.nQuantize > 0 && numel(values) > obj.nQuantize
            mm = minmax(values);
            tresholds = linspace(mm(1), mm(2), obj.nQuantize);
          else
            tresholds = unique(values);
          end
          for treshold = tresholds
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