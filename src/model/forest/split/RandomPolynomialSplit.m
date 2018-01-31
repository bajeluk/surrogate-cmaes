classdef RandomPolynomialSplit < RandomSplit
% RandomPolynomialSplit tries some random polynomial splits and returns the
% best split. It selects a random point from X as origin and generates a
% hyperplane passing through the origin.
  
  properties %(Access = protected)
    split_degree % degree of polynomial
  end
  
  methods
    function obj = RandomPolynomialSplit(options)
      obj = obj@RandomSplit(options);
      obj.split_degree = defopts(options, 'split_degree', 'linear');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [~, d] = size(obj.split_X);
      dPoly = -1;
      for iRepeats = 1:obj.split_nRepeats
        candidate = obj.splitCandidate;
        featuresMin = min(obj.split_X);
        featuresMax = max(obj.split_X);
        switch obj.split_degree
          case 'linear2'
            % each column has values within [min(column), max(column)]
            featuresRand = rand(d, d) ...
              .* repmat(featuresMax - featuresMin, d, 1) ...
              + repmat(featuresMin, d, 1);
            line = featuresRand \ ones(d, 1);
            candidate.splitter = obj.createSplitter(@(X) ...
              X * line - 1);
          otherwise
            % select point where hyperplane passes through
            origin = rand(1, d) ...
              .* (featuresMax - featuresMin) ...
              + featuresMin;
            %origin = datasample(obj.split_X, 1);
            if dPoly < 0
              dPoly = size(generateFeatures(1:d, obj.split_degree, false), 2);
            end
            % select a direction of the hyperplane
            angles = rand(dPoly, 1) * pi - pi/2;
            % convert direction to weights
            weights = tan(angles);
            degree = obj.split_degree;
            candidate.splitter = obj.createSplitter(@(X) ...
              generateFeatures(bsxfun(@minus, X, origin), degree, false)...
              * weights);
        end
        [candidate.gain, candidate.leftID, candidate.rightID] = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end

