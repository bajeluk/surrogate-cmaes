classdef (Abstract) WeakModel < handle
  
  properties
  end
  
  methods (Access = protected)
    function obj = WeakModel(modelOptions)
    end
  end

  methods (Abstract)
    obj = trainModel(obj, X, y)
    % train the model based on the data (X,y)

    [yPred, sd2, ci] = modelPredict(obj, X)
    % predicts the function values in new points X
    
    N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
  end
  
end