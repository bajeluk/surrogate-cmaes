classdef (Abstract) Model
  properties (Abstract)
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration     % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
  end

  methods (Abstract)
    % obj = Model(modelOptions, xMean)
    % constructor

    nData = getNTrainData(obj)
    % returns the required number of data for training the model

    obj = train(obj, X, y, xMean, generation)
    % train the model based on the data (X,y)

    [y, dev] = predict(obj, X)
    % predicts the function values in new points X
  end
end
