classdef GpModel < Model
  properties    % derived from abstract class "Model"
    dim                 % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean           % mean of the generation when the model was built
    dataset             % .X and .y
    useShift = false;
    shiftMean           % vector of the shift in the X-space
    shiftY = 0;         % shift in the f-space
    options

    hyp
    meanFcn
    covFcn
    likFcn
    infFcn
  end

  properties (Access = protected)
    defaultFminconOpts = optimset( ...
      'GradObj', 'on', ...
      'TolFun', 1e-8, ...
      'TolX', 1e-8, ...
      'MaxIter', 1000, ...
      'MaxFunEvals', 1000, ...
      'Display', 'final' ...
      )
  end

  % TODO:
  %   [ ] use the CMA-ES' covariance for Mahalanobis distance

  methods (Access = public)
    function obj = GpModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');
      
      obj.options = modelOptions;
      if (~isempty(modelOptions) && isfield(modelOptions, 'useShift'))
        obj.useShift = modelOptions.useShift;
      end
      obj.dim     = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY  = 0;

      % this is a MOCK/TEST INITIALIZATION!
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
      % TODO
      %   [ ] implement choosing the best covariance function according to
      %       the test ordinal regression capabilities

      assert(size(xMean,1) == 1, 'GpModel.train(): xMean is not a row-vector.');
      obj.trainGeneration = generation;
      obj.trainMean = xMean;
      obj.hyp.mean = mean(y);
      obj.dataset.X = X;
      obj.dataset.y = y;

      if (~isfield(obj.options, 'trainAlgorithm'))
        obj.options.trainAlgorithm = 'minimize';
      end

      if (strcmpi(obj.options.trainAlgorithm, 'minimize')
        % modelTrainNErrors = 0;
        warning('off');
        [hyp_, fX, iters] = minimize(obj.hyp, @gp, -100, @infExactCountErrors, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);
        % DEBUG OUTPUT:
        fprintf('  ... minimize() %f --> %f in %d iterations.\n', fX(1), fX(end), iters);
        warning('on');

        % FIXME: holds for infExact() only -- do not be sticked to infExact!!!
        % nErrors = modelTrainNErrors;
        obj.hyp = hyp_;

      else
        % gp() with linearized version of the hyper-parameters
        f = @(par) linear_gp(par, obj.hyp, @infExactCountErrors, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);

        linear_hyp = unwrap(obj.hyp)';
        l_cov = length(obj.hyp.cov);
        
        % lower and upper bounds
        lb_hyp.cov = -2 * ones(size(obj.hyp.cov));
        lb_hyp.inf = log(1e-7);
        lb_hyp.lik = log(1e-7);
        lb_hyp.mean = -Inf;
        lb = unwrap(lb_hyp)';
        ub_hyp.cov = 25 * ones(size(obj.hyp.cov));
        ub_hyp.inf = log(7);
        ub_hyp.lik = log(7);
        ub_hyp.mean = Inf;
        ub = unwrap(ub_hyp)';

        if (strcmpi(obj.options.trainAlg, 'fmincon'))
         fminconOpts = defaultFminconOpts;
        end
        if (strcmpi(obj.options.trainAlg, 'cmaes'))
        else
          error('GpModel.train(): train algorithm "%s" is not known.\n', obj.options.trainAlgorithm);
        end
      end
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

  methods (Access = private)
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

function [nlZ dnlZ] = linear_gp(linear_hyp, s_hyp, inf, mean, cov, lik, x, y)
  hyp = rewrap(s_hyp, linear_hyp');
  [nlZ s_dnlZ] = gp(hyp, inf, mean, cov, lik, x, y);
  dnlZ = unwrap(s_dnlZ)';
end

function [c ceq] = nonlincons(x)
  % checks if the values x(1:(end-4)) are within 2.5 off median
  MAX_DIFF = 2.5;
  ceq = [];
  assert(size(x,2) == 1, 'Argument for nonlincons is not a vector');
  c = zeros(size(x));
  % test only for covariance parameters
  % TODO: are there always 4 more parameters?!
  c = abs(x(1:end-4) - median(x(1:end-4))) - MAX_DIFF;
end
