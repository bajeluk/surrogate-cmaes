classdef AxisSplit < Split
% AxisSplit finds the best axis parallel split
  
  properties %(Access = protected)
    nQuantize % quantization of tresholds
  end

  methods
    function obj = AxisSplit(options)
      obj = obj@Split(options);
      obj.nQuantize = defopts(options, 'nQuantize', 0);
    end

    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      [n, d] = size(obj.X);
      for feature = 1:d
        featureSelector = (1:d == feature)';
        values = obj.X(:, feature)';
        if obj.nQuantize > 0 && numel(values) > obj.nQuantize
          mm = minmax(values);
          tresholds = linspace(mm(1), mm(2), obj.nQuantize);
        else
          tresholds = unique(values);
        end
        for treshold = tresholds
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
  
end