classdef ResidualObliqueSplit < Split
% ResidualObliqueSplit fits a polynomial model and splits the points into
% two classes based on which side of the line they lie

  properties %(Access = protected)
    modelSpec % degree for fitted polynomial
  end
    
  methods
    function obj = ResidualObliqueSplit(options)
      obj = obj@Split(options);
      obj.modelSpec = defopts(options, 'degree', 'constant');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      trans = obj.transformation;
      candidate = obj.splitCandidate;
      % linear regression
      model = PolynomialModel(struct('modelSpec', obj.modelSpec));
      model = model.trainModel(obj.X, obj.y);
      c = model.modelPredict(obj.X) < obj.y;
      
      switch obj.modelSpec
        case 'constant'
          discrimTypes = {'linear', 'quadratic'};
        case 'linear'
          discrimTypes = {'linear', 'quadratic'};
        case 'quadratic'
          discrimTypes = {'quadratic'};
        otherwise
          discrimTypes = {};
      end
      for i = 1:numel(discrimTypes)
        discrimType = discrimTypes{i};
        try
          model = fitcdiscr(obj.X, c, 'DiscrimType', discrimType);
        catch
          % singular covariance matrix
          pseudoDiscrimType = strcat('pseudo', discrimType);
          model = fitcdiscr(obj.X, c, 'DiscrimType', pseudoDiscrimType);
        end
        candidate.splitter = @(X)...
          model.predict(transformApply(X, trans)) == 1;
        candidate.gain = splitGain.get(candidate.splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end