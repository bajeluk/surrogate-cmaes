classdef RegressPolynomialModel < WeakModel
  
  properties %(Access = protected)
    modelSpec % model specification (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    coeff % coefficients
    coeffCov % coefficient covariance
    features % used features
  end
  
  methods
    function obj = RegressPolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.modelSpec = defopts(modelOptions, 'modelSpec', 'constant');
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      XP = generateFeatures(X, obj.modelSpec, true);
      %warning('off', 'stats:regress:RankDefDesignMat');
      obj.coeff = regress(y, XP);
      %warning('on', 'stats:regress:RankDefDesignMat');
      obj.features = obj.coeff ~= 0;
      obj.coeff = obj.coeff(obj.features);
      XP = XP(:, obj.features);
      M = XP' * XP;
      Mi = inv(M);
      yPred = XP * obj.coeff;
      % var(b) = E(b^2) * (X'*X)^-1
      r = y - yPred;
      mse = r' * r / numel(r);
      obj.coeffCov = mse * Mi;
    end
    
    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      XP = generateFeatures(X, obj.modelSpec, true);
      if ~isempty(obj.features)
        XP = XP(:, obj.features);
      end
      [yPred] = XP * obj.coeff;
      if nargout >= 2
        % sd2 = diag(XP * obj.coeffCov * XP');
        sd2 = sum(XP * obj.coeffCov .* XP, 2);
        if nargout >= 3
          ci = varToConfidence(yPred, sd2);
        end
      end
    end
  end
  
end

