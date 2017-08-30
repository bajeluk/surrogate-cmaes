classdef ConstantModel < WeakModel
  
  properties %(Access = protected)
    coeff % predicted value
    coeffCov % standard deviation
  end
  
  methods
    function obj = ConstantModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.coeff = defopts(modelOptions, 'coeff', NaN);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      if isnan(obj.coeff)
        obj.coeff = mean(y);
        obj.coeffCov = var(y);
      else
        obj.coeffCov = mean((y - obj.coeff).^2);
      end
    end

    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      yPred = repmat(obj.coeff, size(X, 1), 1);
      if nargout >= 2
        sd2 = repmat(obj.coeffCov, size(X, 1), 1);
        if nargout >= 3
          ci = WeakModel.getCi(yPred, sd2);
        end
      end
    end
  end
  
end