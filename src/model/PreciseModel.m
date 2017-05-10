classdef PreciseModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainSigma            % sigma of the generation when the model was built
    trainBD               % BD of the generation when the model was built
    trainMean             % mean of the generation when the model was built
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    options
    predictionType = 'fValues';     % type of prediction (f-values, PoI, EI)
    transformCoordinates = false;   % whether use transformation in the X-space
    stateVariables       % variables needed for sampling new points as CMA-ES do
    sampleOpts           % options and settings for the CMA-ES sampling
    dimReduction          % Reduce dimensionality for model by eigenvectors of covatiance matrix in percentage
    
    bbob_func
  end

  methods
    function obj = PreciseModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');
      
      obj.options   = modelOptions;
      obj.useShift  = defopts(obj.options, 'useShift', false);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      obj.trainBD   = eye(obj.dim);
      obj.dimReduction = 1; %precise model do not use dim reduction

      % BBOB function ID
      % this has to called in opt_s_cmaes due to the speed optimization
      %   handlesF = benchmarks('handles');
      if (isfield(modelOptions, 'bbob_func'))
        obj.bbob_func = modelOptions.bbob_func;
      else
        error('PreciseModel: no BBOB function handle spedified!');
      end
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model

      % BBOB does not need any data, of course, but pretend that it does
      nData = 2 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the GP model based on the data (X,y)

      assert(size(xMean,1) == 1, 'GpModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      % y = (feval(obj.bbob_func, X'))';
      XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
      y = (feval(obj.bbob_func, XWithShift'))';
      y = y + obj.shiftY;
      sd2 = zeros(size(X,1),1);
    end

    function trained = isTrained(obj)
      % the precise model is always trained :)
      trained = true;
    end
  end

end
