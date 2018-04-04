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
      nTreshPerDim = obj.getNTresh();
      for feature = 1:d
        if nTreshPerDim(feature) > 0
          featureSelector = (1:d == feature)';
          values = obj.split_X(:, feature)';
          % calculate tresholds
          tresholds = obj.calcTresholds(values, d, nTreshPerDim(feature));
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
  
  methods (Access = protected)
    function nTreshPerDim = getNTresh(obj)
    % calculates the number of tresholds per dimension according to the 
    % budget of hyperplanes
      [n, dim] = size(obj.split_X);
      % get prescribed number of tresholds
      nTresh = obj.getNQuant(dim);
      % number of tresholds should not be greater than number of points - 1
      nTresh = min(nTresh, n - 1);
      % calculate hyperplane budget
      maxHyp = obj.getMaxHyp(n, dim);
      % too many hyperplanes requires adjustment of treshold numbers
      if dim*nTresh > maxHyp
        nTreshPerDim = floor(maxHyp/dim)*ones(1, dim);
        nRemain = mod(maxHyp, dim);
        % choose dimensions with extra tresholds at random
        nTreshPerDim(randperm(dim)) = nTreshPerDim + [ones(1, nRemain), zeros(1, dim - nRemain)];
      else
        nTreshPerDim = nTresh*ones(1, dim);
      end
    end
  end
  
end