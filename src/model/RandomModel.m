classdef RandomModel < Model
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
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling
    dimReduction          % reduce dimensionality for model by eigenvectors of covariance matrix in percentage
    
    bbob_func
    distribution          % distribution of y-values (uniform, normal, permute)
    yRange                % range of y-values (bbob, dataset)
  end

  methods
    function obj = RandomModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'RandomModel (constructor): xMean is not a row-vector.');
      
      obj.options   = modelOptions;
      obj.useShift  = defopts(obj.options, 'useShift', false);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      % y-values can be shifted (e.g., due to low error stopping criteria)
      obj.shiftY    = defopts(obj.options, 'shiftY', 0);
      obj.trainBD   = eye(obj.dim);
      obj.dimReduction = 1; % random model do not use dim reduction
      % distribution of y-values
      obj.distribution = defopts(obj.options, 'distribution', 'uniform');
      obj.yRange = defopts(obj.options, 'yRange', 'bbob');

      % BBOB function ID
      % this has to called in opt_s_cmaes due to the speed optimization
      %   handlesF = benchmarks('handles');
      if (isfield(modelOptions, 'bbob_func') && strcmp(obj.yRange, 'bbob'))
        obj.bbob_func = modelOptions.bbob_func;
      elseif strcmp(obj.yRange, 'bbob')
        error('RandomModel: no BBOB function handle specified!');
      end
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model

      % BBOB does not need any data, of course, but pretend that it does
      nData = 2 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
      % train the random model based on the data (X,y)

      assert(size(xMean,1) == 1, 'RandomModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.dataset.X = X;
      obj.dataset.y = y;
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      % y = (feval(obj.bbob_func, X'))';
      nPoints = size(X, 1);
      if nPoints < 1
        y = [];
        sd2 = [];
        return
      end
      
      XWithShift = X - repmat(obj.shiftMean, nPoints, 1);
      % get points on range calculation
      if strcmp(obj.yRange, 'bbob')
        yr = (feval(obj.bbob_func, XWithShift'))';
      else
        yr = obj.dataset.y;
      end
      % choose appropriate distribution
      switch obj.distribution
        case {'gauss', 'normal'}
          y = mean(yr) + randn(nPoints, 1) * std(yr);
        case 'uniform'
          y = (max(yr) - min(yr))*rand(nPoints, 1) + min(yr);
        case {'perm', 'permute', 'randperm'}
          if nPoints <= numel(yr)
            y = yr(randperm(numel(yr), nPoints));
          else
            y = yr(randi(numel(yr), nPoints, 1));
          end
      end
      y = y + obj.shiftY;
      sd2 = var(y)*ones(nPoints, 1);
    end

    function trained = isTrained(obj)
      % the random model is always trained :)
      trained = true;
    end
  end

end