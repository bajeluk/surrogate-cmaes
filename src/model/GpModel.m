classdef GpModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration     % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    options

    hyp
    meanFcn
    covFcn
    likFcn
    infFcn
    dataset
  end

  methods
    function obj = GpModel(modelOptions, xMean)
      % constructor
      obj.options = modelOptions;
      obj.dim     = size(xMean, 2);

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

      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.hyp.mean = mean(y);
      obj.dataset.X = X;
      obj.dataset.y = y;

      % modelTrainNErrors = 0;
      hyp = minimize(obj.hyp, @gp, -100, @infExactCountErrors, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);
      % FIXME: holds for infExact() only -- do not be sticked to infExact!!!
      % nErrors = modelTrainNErrors;
      obj.hyp = hyp;
    end

    function [y, dev] = predict(obj, X)
      % predicts the function values in new points X
      [y, dev] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.dataset.X, obj.dataset.y, X);
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
