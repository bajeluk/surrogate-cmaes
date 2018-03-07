classdef PairObliqueSplit < Split
% PairObliqueSplit for each pair of points constructs a normal vector,
% projects all points onto this vector and uses the projected values as
% tresholds for decision to which side the point belongs

  properties %(Access = protected)
    split_nQuantize % quantization of tresholds 
                    %   0 - all tresholds for one pair
                    %   1, 2, 3, ... number of linearly distributed 
                    %   tresholds per pair
    split_pairFcn   % function to compute the number of pairs according to 
                    % the number of points | function handle
  end

  methods
    function obj = PairObliqueSplit(options)
      obj = obj@Split(options);
      obj.split_nQuantize = defopts(options, 'split_nQuantize', 0);
      obj.split_pairFcn   = defopts(options, 'split_pairFcn', @(x) x*log(x) );
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [n, dim] = size(obj.split_X);
      % all combinations of pairs
      pair = nchoosek(1:n, 2);
      nPairs = max(floor(obj.split_pairFcn(n)), 1);
      pairIds = randperm(nchoosek(n, 2), nPairs);
      for p = 1:nPairs
        % normal vector of the hyperplane
        v = obj.split_X(pair(pairIds(p), 1), :) - obj.split_X(pair(pairIds(p), 2), :);
        % project X onto v
        values = (obj.split_X * v')';
        tresholds = obj.calcTresholds(values, dim);
        for treshold = tresholds
          candidate = obj.splitCandidate;
          candidate.splitter = obj.createSplitter(@(X) ...
            X * v' - treshold);
          [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
    end
  end
end