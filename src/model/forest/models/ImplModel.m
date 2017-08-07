classdef (Abstract) ImplModel < handle % Model
  properties
    dim                  % dimension of the input space X (determined from x_mean)
    trainGeneration      % # of the generation when the model was built
    trainMean            % mean of the generation when the model was trained
    trainSigma           % sigma of the generation when the model was trained
    trainBD              % BD of the generation when the model was trained
    dataset              % .X and .y
    useShift             % whether use shift during generationUpdate()
    shiftMean            % vector of the shift in the X-space
    shiftY               % shift in the f-space
    predictionType       % type of prediction (f-values, PoI, EI)
    transformCoordinates % whether use transformation in the X-space
    stateVariables       % variables needed for sampling new points as CMA-ES do
    sampleOpts           % options and settings for the CMA-ES sampling
  end
  
  methods (Access = protected)
    function obj = ImplModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'Model (constructor): xMean is not a row-vector.');
      
      % computed values
      obj.dim       = size(xMean, 2);
      obj.useShift  = defopts(modelOptions, 'useShift', false);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      
      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
    end
  end

  methods (Abstract)
    nData = getNTrainData(obj)
    % returns the required number of data for training the model

    obj = trainModel(obj, X, y, xMean, generation)
    % train the model based on the data (X,y)

    [y, sd2] = modelPredict(obj, X)
    % predicts the function values in new points X
    % returns empty @y on error
  end
end