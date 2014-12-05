classdef (Abstract) Model
  properties (Abstract)
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration     % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    shiftMean           % vector of the shift in the X-space
    shiftY              % shift in the f-space
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
    % returns empty @y on error
  end

  methods  
    function trained = isTrained(obj)
    % check whether the model is already trained
      trained = (obj.trainGeneration >= 0);
    end

    function obj = shift(obj, xMean)
    % transforms the trained model to new coordinates
    % further predictions will be made according to this shift
    % does not re-train the model, but estimates shift in the f-values
    %   according to overlapping regions of the new and old spaces

      % simulate "older" dataset and predict there y-values
      invShiftedDataset = obj.dataset.X - (xMean - obj.trainMean);
      invShiftedY = predict(obj, invShiftedDataset);
      % set the vector of the shift (this has to be after calling
      % predict()!!!)
      obj.shiftMean = xMean - obj.trainMean; 
      % estimate the shift in the new position of the model
      % it is a pair-wise comparison, TODO: is this formula ok?
      obj.shiftY = median(obj.dataset.y - invShiftedY);
    end

    function [obj, counteval, newX, newY] = shiftReevaluate(obj, xMean, xValid, nOrigEvals, fitfun_handle, varargin)
    % transforms the trained model to new coordinates of @xMean
    % further predictions will be made according to this shift
    % does not re-train the model, but estimates the shift in the f-values
    %   according to the evaluation of a point from the Model.dataset (shifted
    %   in the new position)
    % @xValid      -- set of CMA-ES-generated points (from sigma*N(xMean, BD) )
    % @nOrigEvals  -- the number of allowed orig. evaluations for this run, >= 1
    % @varargin    -- arguments for the fitfun_handle function
    % returns:
    %   @counteval -- the number of used original evaluations (non NaN)
    %   @newX      -- coordinates of the new evaluated points
    %   @newY      -- f-values of the new evaluated points
    %
    % FUTURE WORK: 
    % * TODO more points to reevaluate because of conrolling the model precision
    %        (@nOrigEvals parameter)
    % * TODO add a feedback if the model is not sufficiently precise

      % vector of the shift
      obj.shiftMean = xMean - obj.trainMean; 

      y = NaN; x = zeros(1,obj.dim);
      newX = []; newY = [];
      deniedIdxs = [];
      counteval = 0;
      countevalNaN = 0;
      while (isnan(y) && ~isempty(x))
        [x, datasetIdx] = getNearMean(obj, xMean, deniedIdxs);
        if (isempty(x))
          warning(['Model.shiftReevaluate(): fitness in shifted dataset is all NaN. Stopping using model. CountEvalNaN = ' num2str(countevalNaN)]);
          obj = [];
          return;
        end
        invShiftedX = x - obj.shiftMean;
        % evaluate the shifted archive point
        y = feval(fitfun_handle, invShiftedX', varargin{:});
        counteval = counteval + sum(~isnan(y));
        countevalNaN = countevalNaN + sum(isnan(y));
        if (~any(isnan(y)))
          newX = [newX; invShiftedX];
          newY = [newY; y'];
        end
        deniedIdxs = [deniedIdxs datasetIdx];
        % calculate new shiftY
        obj.shiftY = obj.dataset.y(datasetIdx) - y;
        % check that the model is not flat
        yValid = obj.predict(xValid);
        ySort = sort(yValid);
        thirdY = ySort(ceil(1.1+end/3));
        if ((1 - ySort(1)/thirdY) < 1e-2)
          disp('Model.shiftReevaluate(): fitness is flat. Stopping using model.');
          obj = [];
          return;
        end
      end
    end
  end

  methods (Access = protected)
    function [x, datasetIdx] = getNearMean(obj, xMean, deniedIdxs)
    % Returns a point @x from the training dataset (with datasetIdx) near the
    % @trainMean, possibly in the direction of the new @xMean. 
    % This point is to be evaluated and used for shiftReevaluate()
    % @deniedIdxs       vector of indices to the @database which are not
    %                   allowed to use

      % take a point on the line from the trainMean towards the new xMean
      middlePoint = obj.trainMean + 0.25 * (xMean - obj.trainMean);
      % take only not-denied entries in the dataset
      isAllowed = logical(ones(size(obj.dataset.y)));
      isAllowed(deniedIdxs) = 0;
      if (sum(isAllowed) == 0)
        % there is no more allowed dataset entries
        x = []; datasetIdx = [];
        return
      end
      % find the nearest poin to the middlePoint
      sqDists = sum((obj.dataset.X(isAllowed,:) ...
          - repmat(middlePoint, length(obj.dataset.y)-length(deniedIdxs), 1)).^2, 2);
      [~, datasetIdx] = min(sqDists);
      % we got the index only in the "allowed" subset of the dataset
      allowedIdxs = find(isAllowed);

      datasetIdx = allowedIdxs(datasetIdx);
      x = obj.dataset.X(datasetIdx, :);
    end
  end
end
