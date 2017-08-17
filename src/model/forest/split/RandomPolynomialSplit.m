classdef RandomPolynomialSplit < Split
% RandomPolynomialSplit tries some random polynomial splits and returns the
% best split. It selects a random point from X as origin and generates a
% hyperplane passing through the origin.
  
  properties % (Access = protected)
    degree % degree of polynomial
    nRepeats % number of random repeats
  end
  
  methods
    function obj = RandomPolynomialSplit(...
        transformationOptions, degree, nRepeats)
      obj = obj@Split(transformationOptions);
      obj.degree = degree;
      obj.nRepeats = nRepeats;
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
        switch obj.degree
          case 'linear2'
            % each column has values within [min(column), max(column)]
            featuresRand = rand(d, d) ...
              .* repmat(featuresMax - featuresMin, d, 1) ...
              + repmat(featuresMin, d, 1);
            line = featuresRand \ ones(d, 1);
            candidate.splitter = @(X)...
              transformApply(X, trans) * line <= 1;
          otherwise
            % select point where hyperplane passes through
            origin = rand(1, d) ...
              .* (featuresMax - featuresMin) ...
              + featuresMin;
            %origin = datasample(obj.X, 1);
            degree = obj.degree;
            nPoly = size(generateFeatures(1:d, degree, false), 2);
            % select a direction of the hyperplane
            angles = rand(nPoly, 1) * pi - pi/2;
            % convert direction to weights
            weights = tan(angles);
            candidate.splitter = @(X) ...
              generateFeatures(...
                bsxfun(@minus, transformApply(X, trans), origin), ...
                degree, ...
                false ...
              ) * weights <= 0;
        end
        candidate.gain = splitGain.get(splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end

