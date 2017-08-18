classdef RandomRbfSplit < RandomSplit
% RandomRbfSplit tries some random RBF splits and returns the
% best split. It selects a random point from X as origin and generates a
% sphere with given metric

  properties %(Access = protected)
    metric % metric
  end
  
  methods
    function obj = RandomRbfSplit(transformationOptions, nRepeats, ...
        metric)
      obj = obj@RandomSplit(transformationOptions, nRepeats);
      obj.metric = metric;
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for iRepeats = 1:obj.nRepeats
        candidate = obj.splitCandidate;
        featuresMin = min(obj.X);
        featuresMax = max(obj.X);
        % select point where the RBF oiginates
        origin = rand(1, d) ...
          .* (featuresMax - featuresMin) ...
          + featuresMin;
        %origin = datasample(X, 1);
        metric = obj.metric;
        distances = pdist2(X, origin, metric);
        maxDistance = max(distances);
        minDistance = min(distances);
        treshold = rand() * (maxDistance - minDistance) + minDistance;
        switch metric
          case 'mahalanobis'
            C = nancov(X);
            candidate.splitter = @(X)...
              pdist2(transformApply(X, trans), origin, metric, C) ...
              <= treshold;
          otherwise
            candidate.splitter = @(X)...
              pdist2(transformApply(X, trans), origin, metric) ...
              <= treshold;
        end
        candidate.gain = splitGain.get(splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end

