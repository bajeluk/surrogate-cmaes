classdef ResidualObliqueSplit < Split
% ResidualObliqueSplit fits a polynomial model and splits the points into
% two classes based on which side of the line they lie

  properties %(Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    degree % degree for fitted polynomial
  end
    
  methods
    function obj = ResidualObliqueSplit(transformationOptions)
      obj = obj@Split(transformationOptions);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      trans = obj.transformation;
      [n, d] = size(obj.X);
      candidate = obj.splitCandidate;
      % linear regression
      X1 = generateFeatures(obj.X, obj.degree, true);
      w = X1 / obj.y;
      if stricmp(obj.discrimType, obj.degree)
        % classify according to residuals
        candidate.splitter = @(X)...
          generateFeatures(transformApply(X, trans), 'linear', true) ...
          * w < 0;
      else
        % create classes from residuals
        c = X1 * w < obj.y;
        model = fitcdiscr(X, c, 'DiscrimType', obj.discrimType);
        candidate.splitter = @(X)...
          model.predict(transformApply(X, trans)) == 1;
      end
      candidate.gain = splitGain.get(splitter);
      best = candidate;
    end
  end
end