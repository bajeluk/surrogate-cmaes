classdef ForestModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was built
    trainSigma            % sigma of the generation when the model was built
    trainBD               % BD of the generation when the model was built
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling
    options
    
    forestModel           % forest model
  end

  methods
    function obj = ForestModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'ForestModel (constructor): xMean is not a row-vector.');
      obj.options = modelOptions;
      
      % computed values
      obj.useShift  = defopts(obj.options, 'useShift', false);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      
      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
      
      % forest options
      obj.forestModel = defopts(modelOptions, 'forestModel', []);
      if isempty(obj.forestModel)
        modelOptions = struct;
        modelOptions.nTrees = 10;
        obj.forestModel = RandomForestModel(modelOptions);
      end
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: check correctness of the following expression
      nData = 1;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the GP model based on the data (X,y)

      assert(size(xMean,1) == 1, 'ForestModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
      
      obj.forestModel = obj.forestModel.trainModel(X, y);
        
      % count train MSE
      trainMSE = mean((y - obj.forestModel.modelPredict(X)).^2);
      fprintf('  ForestModel: train MSE = %f\n', trainMSE);
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        
        [y, sd2] = obj.forestModel.modelPredict(XWithShift);
        
        y = y + obj.shiftY;
      else
        y = []; sd2 = [];
        fprintf(2, 'RfModel.predict(): the model is not yet trained!\n');
      end
    end
  end
end

