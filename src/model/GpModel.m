classdef GpModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    shiftMean           % vector of the shift in the X-space
    shiftY = 0;         % shift in the f-space
    options

    hyp
    meanFcn
    covFcn
    likFcn
    infFcn
  end

  methods
    function obj = GpModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');
      
      % persistent isInitialized;
      %
      % % intialize GPML library
      % if (isempty(isInitialized) || ~isInitialized)
      %   if (isfield(modelOptions, 'path') && ~isempty(modelOptions.path))
      %     addpath(modelOptions.path);
      %   end
      %   % if (isfield(modelOptions, 'initScript') && ~isempty(modelOptions.initScript))
      %   %   run(modelOptions.initScript);
      %   % end
      %   isInitialized = true;
      % end
      
      obj.options = modelOptions;
      obj.dim     = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY  = 0;

      % this is a MOCK/TEST IMPLEMENTATION!
      obj.hyp.cov = log([0.05 0.1]);
      obj.hyp.inf = log(1e-2);
      obj.hyp.lik = log(0.0001);
      obj.meanFcn = @meanConst;
      obj.covFcn  = @covSEiso;
      obj.likFcn  = @likGauss;
      obj.infFcn  = @infExact;
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = 3 * obj.dim;
    end

    function obj = train(obj, X, y, xMean, generation)
      % train the GP model based on the data (X,y)

      assert(size(xMean,1) == 1, 'GpModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.hyp.mean = mean(y);
      obj.dataset.X = X;
      obj.dataset.y = y;

      % modelTrainNErrors = 0;
      warning('off');
      [hyp_, fX, iters] = minimize(obj.hyp, @gp, -100, @infExactCountErrors, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);
      % DEBUG OUTPUT:
      % fprintf('  minimize() %f --> %f in %d iterations.\n', fX(1), fX(end), iters);
      warning('on');
      % TODO: do not use the model if it is not able to be learn,
      %       i.e. when (fX(1) - fX(end)) is almost zero

      % FIXME: holds for infExact() only -- do not be sticked to infExact!!!
      % nErrors = modelTrainNErrors;
      obj.hyp = hyp_;
    end

    function [y, dev] = predict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        % calculate GP models' prediction in X
        [y, dev] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.dataset.X, obj.dataset.y, XWithShift);
        % apply the shift in the f-space (if there is any)
        y = y + obj.shiftY;
      else
        y = []; dev = [];
        warning('GpModel.predict(): the model is not yet trained!');
      end
    end

    function trained = isTrained(obj)
      % check whether the model is already trained
      trained = (obj.trainGeneration >= 0);
    end
  end

end

function [post, nlZ, dnlZ] = infExactCountErrors(hyp, mean, cov, lik, x, y)
  try
    [post, nlZ, dnlZ] = infExact(hyp, mean, cov, lik, x, y);
  catch err
    % modelTrainNErrors = modelTrainNErrors + 1;
    throw(err);
  end
end
