classdef PolynomialModel < AbstractProbabilityModel
  properties %(Access = protected)
    modelSpec % model specification (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    model % trained model
  end
  
  methods
    function obj = PolynomialModel(modelOptions, xMean)
      % constructor
      obj = obj@AbstractProbabilityModel(modelOptions, xMean);
      % specific model options
      obj.modelSpec = defopts(modelOptions, 'modelSpec', 'constant');
    end
    
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      switch obj.modelSpec
        case 'constant'
          nData = 2;
        case 'linear'
          nData = 3;
        otherwise
          nData = 4;
      end
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      if size(X, 1) == 1
        X = [X; X];
        y = [y; y];
      end
      obj.dataset.X = X;
      obj.dataset.y = y;
      warning('off', 'stats:LinearModel:RankDefDesignMat');
      obj.model = fitlm(X, y, obj.modelSpec);
      warning('on', 'stats:LinearModel:RankDefDesignMat');
    end
    
    
    function [y, sd2, ci, p] = modelPredict(obj, X)
      % predicts the function values in new points X
      [y, ci] = obj.model.predict(X);
      if nargout >= 1
        % sd2 = diag(X * obj.model.CoefficientCovariance * X');
        XP = generateFeatures(X, obj.modelSpec, true, true);
        sd2 = sum(XP * obj.model.CoefficientCovariance .* XP, 2);
        if nargout >= 4
          % assume normal distribution and compute prediction probability
          p = normpdf(y, y, sqrt(sd2));
        end
      end
    end
  end
end

