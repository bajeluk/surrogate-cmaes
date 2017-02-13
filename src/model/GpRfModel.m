classdef GpRfModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling

    % GpRfModel specific fields
    gp                    % gaussian process model
    rf                    % random forest model
    options               % model options

  end

  properties (Access = protected)
  end

  methods (Access = public)
    function obj = GpRfModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpRfModel (constructor): xMean is not a row-vector.');

      % modelOpts structure
      if (isempty(modelOptions))
        obj.options = struct();
      else
        obj.options = modelOptions;
      end
      
      % computed settings
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      obj.useShift  = defopts(obj.options, 'useShift', false);

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      
      % construct models
      obj.rf = RfModel(modelOptions, xMean);
      obj.gp = GpModel(modelOptions, xMean);
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = max(obj.rf.getNTrainData, obj.gp.getNTrainData);
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the GP + RF model based on the data (X,y)

      assert(size(xMean,1) == 1, '  GpRfModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;

      % normalize y if specified
      % TODO: remove?
%       if (obj.options.normalizeY)
%         obj.shiftY = mean(y);
%         obj.stdY  = std(y);
%         yTrain = (y - obj.shiftY) / obj.stdY;
%       else
%         obj.shiftY = 0;
%         obj.stdY  = 1;
%         yTrain = y;
%       end
      
      % train random forest
      obj.rf = obj.rf.trainModel(X, y, xMean, generation);
      % random forest prediction on the train dataset
      rfPred = obj.rf.modelPredict(X);
      % train Gaussian process on residuals
      obj.gp = obj.gp.trainModel(X, y - rfPred, xMean, generation);
      
    end

    function [y, dev] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % calculate RF models' prediction in X
        [yRF, devRF] = obj.rf.modelPredict(X);
        % calculate GP models' prediction in X
        [yGP, devGP] = obj.gp.modelPredict(X);
        % calculate overall prediction
        y = yRF + yGP;
        % TODO: not sure if:
        dev = devRF + devGP;
      else
        y = []; dev = [];
        fprintf(2, 'GpRfModel.predict(): the model is not yet trained!\n');
      end
    end

  end
end
