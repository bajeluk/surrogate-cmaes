classdef LmfitPolynomialModel < WeakModel
  
  properties %(Access = protected)
    weak_modelSpec % model specification (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    weak_model     % trained model
  end
  
  methods
    function obj = LmfitPolynomialModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      % specific model options
      obj.weak_modelSpec = defopts(modelOptions, 'weak_modelSpec', 'constant');
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)
      if size(X, 1) == 1
        X = [X; X];
        y = [y; y];
      end
      warning('off', 'stats:LinearModel:RankDefDesignMat');
      obj.weak_model = fitlm(X, y, obj.weak_modelSpec);
      warning('on', 'stats:LinearModel:RankDefDesignMat');
    end
    
    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      [yPred, ci] = obj.weak_model.predict(X);
      if nargout >= 2
        sd2 = confidenceToVar(ci);
      end
    end
  end
  
end

