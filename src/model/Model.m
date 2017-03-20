classdef (Abstract) Model
  properties (Abstract)
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

  methods (Abstract)
    % obj = Model(modelOptions, xMean)
    % constructor

    nData = getNTrainData(obj)
    % returns the required number of data for training the model

    obj = trainModel(obj, X, y, xMean, generation)
    % train the model based on the data (X,y)

    [y, sd2] = modelPredict(obj, X)
    % predicts the function values in new points X
    % returns empty @y on error
  end

  methods  
    function trained = isTrained(obj)
    % check whether the model is already trained
      trained = (obj.trainGeneration >= 0);
    end

    function n = getTrainsetSize(obj)
      n = size(obj.getDataset_y(), 1);
    end

    function y = getDataset_X(obj)
      y = obj.dataset.X;
    end

    function y = getDataset_y(obj)
      y = obj.dataset.y;
    end

    function obj = shift(obj, xMean)
    % transforms the trained model to new coordinates
    % further predictions will be made according to this shift
    % does not re-train the model, but estimates shift in the f-values
    %   according to overlapping regions of the new and old spaces

      % simulate "older" dataset and predict there y-values
      invShiftedDataset = obj.getDataset_X() - (xMean - obj.trainMean);
      obj.shiftMean = zeros(size(xMean));
      invShiftedY = predict(obj, invShiftedDataset);
      % set the vector of the shift (this has to be after calling
      % predict()!!!)
      obj.shiftMean = xMean - obj.trainMean; 
      % estimate the shift in the new position of the model
      % it is a pair-wise comparison, TODO: is this formula ok?
      obj.shiftY = median(obj.getDataset_y() - invShiftedY);
    end

    function [obj, counteval, newX, newY, newZ] = generationUpdate(obj, xMean, xValid, zValid, nOrigEvals, fitfun_handle, varargin)
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
      if (obj.useShift)
        obj.shiftMean = xMean - obj.trainMean;
      else
        obj.shiftMean = zeros(size(xMean));
      end

      y = NaN; x = zeros(1,obj.dim);
      newX = []; newY = []; newZ = [];
      counteval = 0;
      countevalNaN = 0;

      if (isempty(nOrigEvals) || (nOrigEvals <= 0))
        % we are not allowed to use any original evaluations
        return;
      end

      % !!!! TODO !!!!
      % !!!! add the newly evaluated points to the ARCHIVE !!!!
      % !!!! TODO !!!!

      deniedIdxs = [];
      while (isnan(y) && ~isempty(x))
        [x, datasetIdx] = getNearMean(obj, xMean, deniedIdxs);
        if (isempty(x))
          warning(['Model.generationUpdate(): fitness in shifted dataset is all NaN. Stopping using model. CountEvalNaN = ' num2str(countevalNaN)]);
          obj = [];
          return;
        end
        % possibly shift 'x' in the X-space
        % 'shiftMean' can be 0 (if obj.useShift == false)
        shiftedX = x + obj.shiftMean;
        % evaluate the (shifted) archive point
        y = feval(fitfun_handle, shiftedX', varargin{:});
        counteval = counteval + sum(~isnan(y));
        countevalNaN = countevalNaN + sum(isnan(y));
        if (all(~isnan(y)))
          newX = [newX; shiftedX];
          newY = [newY; y'];
          newZ = [newZ; zeros(1,obj.dim)];
        end
        deniedIdxs = [deniedIdxs datasetIdx];
        % calculate new shiftY
        dataset_y = obj.getDataset_y();
        obj.shiftY = dataset_y(datasetIdx) - y;
        % check that the model is not flat
        yValid = obj.predict(xValid);
        ySort = sort(yValid);
        thirdY = ySort(ceil(1.1+end/3));
        if ((thirdY - ySort(1)) < 1e-8)
          % disp('Model.generationUpdate(): fitness is flat. Stopping using model.');
          obj = [];
          return;
        end
      end
      % Evaluate other validating points, if there are any
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
          % disp('Model.generationUpdate(): The model is not precise enough, skipping using the model.');
          obj = [];
          return;
        end
      end
    end

    function [y, sd2] = predict(obj, X)
    % Predicts the function values in new points X.

      % transform input variables using Mahalanobis distance
      if obj.transformCoordinates
        % compute coordinates in the (sigma*BD)-basis
        % BDinv = inv(sigma*BD);
        XTransf = ( (obj.trainSigma * obj.trainBD) \ X')';
      else
        XTransf = X;
      end

      % dimensionality reduction
      if (isprop(obj, 'dimReduction') && (obj.dimReduction ~= 1))
        cntDimension = ceil(obj.dim * obj.dimReduction);
        obj.shiftMean = obj.shiftMean(1:cntDimension);
        XtransfReduce = obj.reductionMatrix * XTransf';
        XtransfReduce = XtransfReduce';
      else
        XtransfReduce = XTransf;
      end

      [y, sd2] = modelPredict(obj,XtransfReduce);

    end

    function [output, y] = getModelOutput(obj,X)
    % Predicts the function values, the variance, the probability of 
    % improvement, the expected improvement or the least confidence bound
    % in new points X. Values of PoI  and EI can be transformed to last 
    % known fvalue range.

      [y, sd2] = obj.predict(X);
      fmin = min(obj.getDataset_y());
      fmax = max(obj.getDataset_y());

      switch lower(obj.predictionType)
        case 'fvalues' % mean function values
          output = y;

        case 'sd2' % variance
          output = sd2;

        case 'poi' % probability of improvement
          target = fmin - 0.05 * (fmax - fmin);
          output = getPOI(X, y, sd2, target);

        case 'ei' % expected improvement
          output = getEI(X, y, sd2, fmin);

        case 'lcb' % lower confidence bound
          output = y - 2 * sqrt(sd2);

        case 'fpoi' % PoI scaled using function values
          target = fmin - 0.05 * (fmax - fmin);
          poi = getPOI(X, y, sd2, target);
          poiMax = max(poi);
          poiMin = min(poi);
          % map the higest PoI to the smallest function value and vice versa
          output = (fmax-fmin)*(poiMax - poi)/(poiMax-poiMin)+fmin;

        case 'fei' % EI scaled using function values
          ei = getEI(X, y, sd2, fmin);
          eiMax = max(ei);
          eiMin = min(ei);
          % map the higest EI to the smallest function value and vice versa
          output = (fmax-fmin)*(eiMax-ei)/(eiMax-eiMin)+fmin;

        otherwise % otherwise return sd2
          output = sd2;
      end

    end

    function obj = train(obj, X, y, stateVariables, sampleOpts, archive, population)
    % train the model based on the data (X,y)
    % if archive is passed and trainsetType is not 'parameters',
    % X and y are ignored and new values for X, y are retrieved from archive 
    % according to the model settings

      xMean = stateVariables.xmean';
      generation = stateVariables.countiter;
      sigma = stateVariables.sigma;
      lambda = stateVariables.lambda;
      BD = stateVariables.BD;
      obj.sampleOpts = sampleOpts;
      obj.trainMean = xMean;
      obj.stateVariables = stateVariables;

      trainsetType = defopts(obj.options, 'trainsetType', 'parameters');
      obj.options.trainsetSizeMax = defopts(obj.options, 'trainsetSizeMax', 15*obj.dim);
      obj.options.trainRange = defopts(obj.options, 'trainRange', 1);
      if (~strcmpi(trainsetType, 'parameters') && exist('archive','var'))
        [X, y] = obj.generateDataset(archive, population);
      end

      if (isempty(X))
        obj.trainGeneration= - 1;
        warning('Model.train() - empty trainset. Considering the model as untrained.');
      else
        % minimal difference between minimal and maximal returned
        % value to regard the model as trained; otherwise, the
        % constant response is mark of a badly trained model
        % and therefor it is marked as untrained
        MIN_RESPONSE_DIFFERENCE = min(1e-8, 0.05 * (max(y) - min(y)));

        % transform input variables using Mahalanobis distance
        if obj.transformCoordinates
          % compute coordinates in the (sigma*BD)-basis
          obj.trainSigma = sigma;
          obj.trainBD = BD;
          XTransf = ( (sigma * BD) \ X')';
        else
          XTransf = X;
        end

        % dimensionality reduction
        if (isprop(obj, 'dimReduction') && (obj.dimReduction ~= 1))
          cntDimension = ceil(obj.dim * obj.dimReduction);
          obj.shiftMean = obj.shiftMean(1:cntDimension);
          changeMatrix = (eye(obj.dim) / BD);
          changeMatrix = changeMatrix(1:cntDimension,:);
          obj.reductionMatrix = changeMatrix;
          XtransfReduce = changeMatrix * XTransf';
          XtransfReduce = XtransfReduce';
        else
          XtransfReduce=XTransf;
        end

        obj = trainModel(obj, XtransfReduce, y, xMean, generation);

        if (obj.isTrained())
          % Test that we don't have a constant model
          [~, xTestValid] = sampleCmaesNoFitness(sigma, lambda, stateVariables, sampleOpts);
          yPredict = obj.predict(xTestValid');
          if (max(yPredict) - min(yPredict) < MIN_RESPONSE_DIFFERENCE)
            fprintf('Model.train(): model output is constant (diff=%e), considering the model as un-trained.\n', max(yPredict) - min(yPredict));
            obj.trainGeneration = -1;
          end
        end
      end
    end
  end

  methods (Access = protected)
    function [x, datasetIdx] = getNearMean(obj, xMean, deniedIdxs)
    % Returns a point @x from the training dataset (with datasetIdx) near the
    % @trainMean, possibly in the direction of the new @xMean. 
    % This point is to be evaluated and used for generationUpdate()
    % @deniedIdxs       vector of indices to the @database which are not
    %                   allowed to use

      % take a point on the line from the trainMean towards the new xMean
      middlePoint = obj.trainMean + 0.25 * (xMean - obj.trainMean);
      % take only not-denied entries in the dataset
      isAllowed = logical(true(size(obj.getDataset_y())));
      isAllowed(deniedIdxs) = 0;
      if (sum(isAllowed) == 0)
        % there is no more allowed dataset entries
        x = []; datasetIdx = [];
        return
      end
      % find the nearest poin to the middlePoint
      dataset_X = obj.getDataset_X();
      sqDists = sum((dataset_X(isAllowed,:) ...
          - repmat(middlePoint, length(obj.getDataset_y())-length(deniedIdxs), 1)).^2, 2);
      [~, datasetIdx] = min(sqDists);
      % we got the index only in the "allowed" subset of the dataset
      allowedIdxs = find(isAllowed);

      datasetIdx = allowedIdxs(datasetIdx);
      x = dataset_X(datasetIdx, :);
    end


    function [X,y] = generateDataset(obj, archive, population)
      xMean = obj.stateVariables.xmean';
      sigma = obj.stateVariables.sigma;
      BD = obj.stateVariables.BD;
      dim = obj.dim;
      [X,y] = archive.getTrainsetData(obj.options.trainsetType,...
          myeval(obj.options.trainsetSizeMax), xMean, obj.options.trainRange,...
          sigma, BD, population);
    end
  end
end
