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

    function [obj, counteval, newX, newY, newZ] = shiftReevaluate(obj, xMean, xValid, zValid, nOrigEvals, fitfun_handle, varargin)
    % transforms the trained model to new coordinates of @xMean
    % further predictions will be made according to this shift
    % does not re-train the model, but estimates the shift in the f-values
    %   according to the evaluation of a point from the Model.dataset (shifted
    %   in the new position)
    % @xValid      -- set of CMA-ES-generated points (from sigma*N(xMean, BD) )
    % @zValid      -- set of CMA-ES-generated points (from N(0,1) )
    % @nOrigEvals  -- the number of allowed orig. evaluations for this run, >= 1
    % @varargin    -- arguments for the fitfun_handle function
    % returns:
    %   @obj       -- is empty if the model is not sufficient for prediction
    %   @counteval -- the number of used original evaluations (non NaN)
    %   @newX      -- new evaluated points (starting from the second, 
    %                 they should be copied from xValid)
    %   @newY      -- f-values of the new evaluated points
    %   @newZ      -- correspoing vectors from normal N(0,1), the first row, which
    %                 corresponds to the shifted archive point, will not be computed!
    %
    % FUTURE WORK: 
    % * TODO more points to reevaluate because of conrolling the model precision
    %        (@nOrigEvals parameter)
    % * TODO add a feedback if the model is not sufficiently precise

      % vector of the shift
      obj.shiftMean = xMean - obj.trainMean; 

      y = NaN; x = zeros(1,obj.dim);
      newX = []; newY = []; newZ = [];
      counteval = 0;
      countevalNaN = 0;
      deniedIdxs = [];
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
        if (all(~isnan(y)))
          newX = [newX; invShiftedX];
          newY = [newY; y'];
          newZ = [newZ; zeros(1,obj.dim)];
        end
        deniedIdxs = [deniedIdxs datasetIdx];
        % calculate new shiftY
        obj.shiftY = obj.dataset.y(datasetIdx) - y;
        % check that the model is not flat
        yValid = obj.predict(xValid);
        ySort = sort(yValid);
        thirdY = ySort(ceil(1.1+end/3));
        if ((thirdY - ySort(1)) < 1e-8)
          % disp('Model.shiftReevaluate(): fitness is flat. Stopping using model.');
          obj = [];
          return;
        end
      end
      % Evaluate validating points
      nOrigEvals = min((nOrigEvals - 1), size(xValid, 1));
      if (nOrigEvals > 0)
        xValid2 = xValid(1:nOrigEvals,:) - repmat(obj.shiftMean,nOrigEvals,1);
        for i = 1:nOrigEvals
          % evaluate the shifted archive point
          y = feval(fitfun_handle, xValid2(i,:)', varargin{:});
          counteval = counteval + sum(~isnan(y));
          countevalNaN = countevalNaN + sum(isnan(y));
          if (all(~isnan(y)))
            newX = [newX; xValid2(i,:)];
            newY = [newY; y'];
            newZ = [newZ; zValid(i,:)];
          end
        end
        % calculate how good the model estimates ordering on the validating set
        yPredict = obj.predict(newX);
        kendall = corr(yPredict, newY, 'type', 'Kendall');
        % TODO: this test is a rule of thumb!!! Test it!
        if (~isnan(kendall)  &&  kendall < 0.3)
          % disp('Model.shiftReevaluate(): The model is not precise enough, skipping using the model.');
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
