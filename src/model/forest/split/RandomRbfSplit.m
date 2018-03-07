classdef RandomRbfSplit < RandomSplit
% RandomRbfSplit tries some random RBF splits and returns the
% best split. It selects a random point from X as origin and generates a
% sphere with given metric.

  properties %(Access = protected)
    split_randrbf_metric % metric (according to the matlab pdist2 
                         % distances)
  end
  
  methods
    function obj = RandomRbfSplit(options)
      obj = obj@RandomSplit(options);
      obj.split_randrbf_metric = defopts(options, 'split_randrbf_metric', 'euclidean');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [~, d] = size(obj.split_X);
      for iRepeats = 1:obj.split_nRepeats
        candidate = obj.splitCandidate;
        featuresMin = min(obj.split_X);
        featuresMax = max(obj.split_X);
        % select point where the RBF originates
        origin = rand(1, d) ...
          .* (featuresMax - featuresMin) ...
          + featuresMin;
        % origin = datasample(X, 1);
        metric = obj.split_randrbf_metric;
        distances = pdist2(obj.split_X, origin, metric);
        maxDistance = max(distances);
        minDistance = min(distances);
        treshold = rand() * (maxDistance - minDistance) + minDistance;
        switch metric
          case 'mahalanobis'
            C = nancov(X);
            candidate.splitter = obj.createSplitter(@(X) ...
              pdist2(X, origin, metric, C) - treshold);
          otherwise
            candidate.splitter = obj.createSplitter(@(X) ...
              pdist2(X, origin, metric) - treshold);
        end
        [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end

