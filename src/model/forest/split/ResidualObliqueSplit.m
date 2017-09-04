classdef ResidualObliqueSplit < Split
% ResidualObliqueSplit fits a polynomial model and splits the points into
% two classes based on which side of the line they lie

  properties %(Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    degree % degree for fitted polynomial
  end
    
  methods
    function obj = ResidualObliqueSplit(options)
      obj = obj@Split(options);
      obj.degree = defopts(options, 'degree', 'linear');
      obj.discrimType = defopts(options, 'discrimType', 'linear');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      trans = obj.transformation;
      [n, d] = size(obj.X);
      candidate = obj.splitCandidate;
      % linear regression
      X1 = generateFeatures(obj.X, obj.degree, true);
      w = X1 / obj.y;
      % create classes from residuals
      c = X1 * w < obj.y;
      try
        model = fitcdiscr(obj.X, c, 'DiscrimType', obj.discrimType);
      catch
        % singular covariance matrix
        pseudoDiscrimType = strcat('pseudo', obj.discrimType);
        model = fitcdiscr(obj.X, c, 'DiscrimType', pseudoDiscrimType);
      end
      candidate.splitter = @(X)...
        model.predict(transformApply(X, trans)) == 1;
      candidate.gain = splitGain.get(splitter);
      best = candidate;
    end
  end
end