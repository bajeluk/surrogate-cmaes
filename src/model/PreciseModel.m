classdef PreciseModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    shiftMean           % vector of the shift in the X-space
    shiftY = 0;         % shift in the f-space
    options
    
    bbob_func
  end

  methods
    function obj = PreciseModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');
      
      obj.options = modelOptions;
      obj.dim     = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY  = 0;

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

    function obj = train(obj, X, y, xMean, generation)
      % train the GP model based on the data (X,y)

      assert(size(xMean,1) == 1, 'GpModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
    end

    function [y, dev] = predict(obj, X)
      % predicts the function values in new points X
      y = (feval(obj.bbob_func, X'))';
      dev = zeros(size(X,1),1);
    end

    function trained = isTrained(obj)
      % the precise model is always trained :)
      trained = true;
    end
  end

end
