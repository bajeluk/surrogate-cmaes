classdef PolynomialModel < ImplModel
  properties %(Access = protected)
    modelSpec % model specification (https://www.mathworks.com/help/stats/fitlm.html#inputarg_modelspec)
    model % trained model
  end
  
  methods
    function obj = PolynomialModel(modelOptions, xMean)
      % constructor
      obj = obj@ImplModel(modelOptions, xMean);
      % specific model options
      obj.modelSpec = defopts(modelOptions, 'modelSpec', 'constant');
    end
    
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = 1;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      warning('off', 'stats:LinearModel:RankDefDesignMat');
      obj.model = fitlm(X, y, obj.modelSpec);
      warning('on', 'stats:LinearModel:RankDefDesignMat');
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      [y, ci] = obj.model.predict(X);
      sd2 = (ci(:,2) - ci(:,1)).^2; % TODO: can we get SD2 from confidence interval?
      sd2 = repmat(obj.model.RMSE^2, size(y, 1), 1);
    end
  end
end

