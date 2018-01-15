classdef RegressPolynomialModel < WeakModel
  
  properties %(Access = protected)
    weak_modelSpec % model specification from MATLAB fitlm function
    % (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    % except 'polyijk' settings
    weak_coeff % coefficients
    weak_coeffCov % coefficient covariance
    weak_features % used features
  end
  
  methods
    function obj = RegressPolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_modelSpec = defopts(modelOptions, 'weak_modelSpec', 'constant');
      modelSpec_types = {'constant', 'linear', 'interactions', 'purequadratic', 'quadratic'};
      assert(any(strcmp(obj.weak_modelSpec, modelSpec_types)), ...
        'Model ''%s'' cannot be specified as weak_modelSpec property for RegressPolynomialModel', ...
        obj.weak_modelSpec);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      XP = generateFeatures(X, obj.weak_modelSpec, true);
      warning('off', 'stats:regress:RankDefDesignMat');
      obj.weak_coeff = regress(y, XP);
      warning('on', 'stats:regress:RankDefDesignMat');
      obj.weak_features = obj.weak_coeff ~= 0;
      obj.weak_coeff = obj.weak_coeff(obj.weak_features);
      XP = XP(:, obj.weak_features);
      M = XP' * XP;
      Mi = inv(M);
      yPred = XP * obj.weak_coeff;
      % var(b) = E(b^2) * (X'*X)^-1
      r = y - yPred;
      mse = r' * r / numel(r);
      obj.weak_coeffCov =  mse * Mi;
    end
    
    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      XP = generateFeatures(X, obj.weak_modelSpec, true);
      if ~isempty(obj.weak_features)
        XP = XP(:, obj.weak_features);
      end
      [yPred] = XP * obj.weak_coeff;
      if nargout >= 2
        % sd2 = diag(XP * obj.weak_coeffCov * XP');
        sd2 = sum(XP * obj.weak_coeffCov .* XP, 2);
        if nargout >= 3
          ci = varToConfidence(yPred, sd2);
        end
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      switch obj.weak_modelSpec
        case 'constant'
          N = ones(size(dim));
        case 'linear'
          N = 1 + dim;
        case 'interactions'
          N = 1 + dim + dim.*(dim-1)/2;
        case 'purequadratic'
          N = 1 + 2*dim;
        case 'quadratic'
          N = 1 + 2*dim + dim.*(dim-1)/2;
        otherwise
          N = [];
      end
    end
    
  end
  
end

