classdef ConstantModel < WeakModel
  
  properties %(Access = protected)
    weak_coeff    % predicted value
    weak_coeffCov % standard deviation
  end
  
  methods
    function obj = ConstantModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_coeff = defopts(modelOptions, 'weak_coeff', NaN);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      if isnan(obj.weak_coeff)
        obj.weak_coeff = sum(y) / numel(y);
      end
      r = y - obj.weak_coeff;
      obj.weak_coeffCov = r' * r / numel(r);
    end

    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      yPred = repmat(obj.weak_coeff, size(X, 1), 1);
      if nargout >= 2
        sd2 = repmat(obj.weak_coeffCov, size(X, 1), 1);
        if nargout >= 3
          ci = varToConfidence(yPred, sd2);
        end
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      N = ones(size(dim));
    end
  end
  
end