classdef AxisSplit < Split
% AxisSplit finds the best axis parallel split
  
  properties %(Access = protected)
    split_nQuantize % quantization of tresholds 
                    %   0 - all tresholds for one dimension
                    %   1, 2, 3, ... number of linearly distributed 
                    %   tresholds per dimension
  end

  methods
    function obj = AxisSplit(options)
      obj = obj@Split(options);
      obj.split_nQuantize = defopts(options, 'split_nQuantize', 0);
    end

    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [~, d] = size(obj.split_X);
      for feature = 1:d
        featureSelector = (1:d == feature)';
        values = obj.split_X(:, feature)';
        tresholds = obj.calcTresholds(values, d);
        % calculate gain for each treshold
        for treshold = tresholds
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
  
end