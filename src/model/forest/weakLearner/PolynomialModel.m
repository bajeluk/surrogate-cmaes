classdef PolynomialModel < WeakModel
  
  properties %(Access = protected)
    weak_modelSpec % model specification (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    weak_coeff % coefficients
    weak_coeffCov % coefficient covariance
    weak_features % used features
  end
  
  methods
    function obj = PolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_modelSpec = defopts(modelOptions, 'weak_modelSpec', 'constant');
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      XP = generateFeatures(X, obj.weak_modelSpec, true);
      M = XP' * XP;
      % check rank deficiency
      [U, S, V] = svd(M, 'econ');
      s = diag(S);
      tol = max(size(M)) * eps(max(s));
      rank = sum(s > tol);
      if rank < size(M, 2)
        % remove dependent columns
        % [~, obj.weak_features] = rref(M);
        [~, R, E] = qr(M, 0);
        s = abs(diag(R));
        tol = max(size(M)) * eps(max(s));
        obj.weak_features = E(s > tol);
        XP = XP(:, obj.weak_features);
        M = M(obj.weak_features, obj.weak_features);
        [U, S, V] = svd(M, 'econ');
      end
      %warning('off', 'MATLAB:rankDeficientMatrix');
      %warning('off', 'MATLAB:singularMatrix');
      warning('off', 'MATLAB:nearlySingularMatrix');
      %Mi = inv(M);
      Mi = V / S * U';
      %obj.weak_coeff = Mi * XP' * y;
      %obj.weak_coeff = M \ (XP' * y);
      obj.weak_coeff = XP \ y;
      %warning('on', 'MATLAB:rankDeficientMatrix');
      %warning('on', 'MATLAB:singularMatrix');
      warning('on', 'MATLAB:nearlySingularMatrix');
      yPred = XP * obj.weak_coeff;
      % var(b) = E(b^2) * (X'*X)^-1
      r = y - yPred;
      mse = r' * r / numel(r);
      obj.weak_coeffCov = mse * Mi;
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
  end
  
end

