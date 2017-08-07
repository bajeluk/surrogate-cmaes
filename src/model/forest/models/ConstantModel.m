classdef ConstantModel < ImplModel
  properties %(Access = protected)
    y % predicted value
    sd2 % standard deviation
  end
  
  methods
    function obj = ConstantModel(modelOptions, xMean)
      % constructor
      obj = obj@ImplModel(modelOptions, xMean);
      % specific model options
      obj.y = defopts(modelOptions, 'y', NaN);
    end
    
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = 0;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the model based on the data (X,y)
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      if isnan(obj.y)
        obj.y = mean(y);
      end
      obj.sd2 = sum((y - obj.y).^2) / size(y, 1);
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      y = repmat(obj.y, size(X, 1), 1);
      sd2 = repmat(obj.sd2, size(X, 1), 1);
    end
  end
  
end