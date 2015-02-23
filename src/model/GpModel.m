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
    nErrors
    trainLikelihood
  end

  properties (Access = protected)
  end

  % TODO:
  %   [ ] use the CMA-ES' covariance for Mahalanobis distance

  methods (Access = public)
    function obj = GpModel(modelOptions, xMean)
      % constructor
      assert(size(xMean,1) == 1, 'GpModel (constructor): xMean is not a row-vector.');

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
      obj.trainLikelihood = Inf;
      obj.useShift  = defopts(obj.options, 'useShift', false);

      % Optimization Toolbox check
      obj.options.trainAlgorithm = defopts(obj.options, 'trainAlgorithm', 'fmincon');
      if (strcmpi(obj.options.trainAlgorithm, 'fmincon') ...
          && ~license('checkout', 'optimization_toolbox'))
        warning('GpModel: Optimization Toolbox license not available. Switching to minimize().');
        obj.options.trainAlgorithm = 'minimize';
      end

      % GP hyper-parameter settings
      if (~isfield(obj.options, 'hyp'))
        obj.options.hyp = struct();
      end
      obj.hyp.inf = defopts(obj.options.hyp, 'inf', log(1e-2)); % should be roughly somewhere between log(1e-3) and log(1e-2)
      obj.hyp.lik = defopts(obj.options.hyp, 'lik', log(0.01));  % should be somewhere between log(0.01) and log(1)
      obj.hyp.cov = defopts(obj.options.hyp, 'cov', log([0.5; 2]));   % should be somewhere between log([0.1 2]) and log([2 1e6])
      covFcn = defopts(obj.options, 'covFcn',  '{@covMaterniso, 5}');
      if (exist(covFcn) == 2)
        % string with name of an m-file function
        obj.covFcn  = str2func(covFcn);
      else
        % something more complex, like '{@covMaterniso, 3}'
        obj.covFcn  = eval(covFcn);
      end
      obj.meanFcn = str2func(defopts(obj.options, 'meanFcn', 'meanConst'));
      obj.likFcn  = str2func(defopts(obj.options, 'likFcn',  'likGauss'));
      obj.infFcn  = str2func(defopts(obj.options, 'infFcn',  'infExactCountErrors'));

      % use POI or EI?
      obj.options.usePOI = defopts(obj.options, 'usePOI', false);
      obj.options.useEI  = defopts(obj.options, 'useEI',  false);
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
      global modelTrainNErrors;

      assert(size(xMean,1) == 1, '  GpModel.train(): xMean is not a row-vector.');
      obj.trainMean = xMean;
      obj.hyp.mean = mean(y);
      obj.dataset.X = X;
      obj.dataset.y = y;

      if (~isfield(obj.options, 'trainAlgorithm'))
        obj.options.trainAlgorithm = 'minimize';
      end
      alg = obj.options.trainAlgorithm;

      if (strcmpi(alg, 'minimize'))
        [obj, fval] = obj.trainMinimize(X, y);
        if (fval < Inf)
          obj.trainGeneration = generation;
        else
          obj.trainGeneration = -1;
        end

      elseif (strcmpi(alg, 'fmincon') ...
              || strcmp(alg, 'cmaes'))
        % gp() with linearized version of the hyper-parameters
        f = @(par) linear_gp(par, obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);

        linear_hyp = unwrap(obj.hyp)';
        l_cov = length(obj.hyp.cov);

        % lower and upper bounds
        [lb_hyp, ub_hyp] = obj.defaultLUBounds();
        lb = unwrap(lb_hyp)';
        ub = unwrap(ub_hyp)';
        opt = [];

        if (strcmpi(alg, 'fmincon'))
          [obj, opt, trainErr] = obj.trainFmincon(linear_hyp, X, y, lb, ub, f);

          if (trainErr)
            disp('Trying CMA-ES...');
            alg = 'cmaes';
          end
        end
        if (strcmpi(alg, 'cmaes'))
          [obj, opt, trainErr] = obj.trainCmaes(linear_hyp, X, y, lb, ub, f);
          if (trainErr)
            return;
          end
        end

        obj.trainGeneration = generation;
        obj.hyp = rewrap(obj.hyp, opt);

        % DEBUG OUTPUT:
        fprintf('.. model-training likelihood = %f\n', obj.trainLikelihood);
        % disp(obj.hyp);
      else
        error('GpModel.train(): train algorithm "%s" is not known.\n', alg);
      end
    end

    function [y, dev] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        % calculate GP models' prediction in X
        [y, dev] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.dataset.X, obj.dataset.y, XWithShift);
        % apply the shift in the f-space (if there is any)
        y = y + obj.shiftY;

        % % Calculate POI if it should be used
        % if (obj.options.usePOI)
        %   % return -POI , because the smaller y-value (bigger probability) the better
        %   y = - getPOI(X, y, dev, min(obj.dataset.y));
        %   dev = zeros(size(dev));
        % end
        %
        % % Calculate EI if it should be used
        % if (obj.options.useEI)
        %   % EI should be negative in promising regions, the lower the better
        %   y = getEI(X, y, dev, min(obj.dataset.y));
        %   dev = zeros(size(dev));
        % end
      else
        y = []; dev = [];
        fprintf(2, 'GpModel.predict(): the model is not yet trained!\n');
      end
    end

    function trained = isTrained(obj)
      % check whether the model is already trained
      trained = (obj.trainGeneration >= 0);
    end
  end

  methods (Access = private)
    function [obj, fval] = trainMinimize(obj, X, y)
      % train the GP model using Rasmussen's minimize() function
      %
      global modelTrainNErrors;

      modelTrainNErrors = 0;
      fval = Inf;
      warning('off');
      try
        [hyp_, fval, iters] = minimize(obj.hyp, @gp, -100, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, X, y);
      catch
        fprintf(2, 'minimize() failed.\n');
        warning('on');
        obj.nErrors = modelTrainNErrors;
        return;
      end
      % DEBUG OUTPUT:
      fprintf('  ... minimize() %f --> %f in %d iterations.\n', fval(1), fval(end), iters);
      warning('on');

      obj.nErrors = modelTrainNErrors;
      obj.trainLikelihood = fval(end);
      obj.hyp = hyp_;
    end

    function [obj, opt, trainErr] = trainFmincon(obj, linear_hyp, X, y, lb, ub, f);
      % train with Matlab's fmincon() from the Optimization toolbox
      %
      global modelTrainNErrors;
      trainErr = false;
      opt = [];

      [fminconOpts, nonlnc] = obj.defaultFminconOpts();
      try
        initial = f(linear_hyp');
      catch err
        initial = NaN;
      end
      if isnan(initial)
        % the initial point is not valid
        disp('  GpModel.train(): fmincon -- initial point is not valid.');
        trainErr = true;
      else
        % training itself
        disp(['Model training (fmincon), init fval = ' num2str(initial)]);
        try
          modelTrainNErrors = 0;
          [opt, fval] = fmincon(f, linear_hyp', [], [], [], [], lb, ub, nonlnc, fminconOpts);
          obj.nErrors = modelTrainNErrors;
          obj.trainLikelihood = fval;
          if (isnan(fval))
            trainErr = true;
          end
        catch err
          obj.nErrors = modelTrainNErrors;
          fprintf(2, '  GpModel.train() ERROR: fmincon() ended with an exception: %s\n', err.message);
          trainErr = true;
        end
      end
    end

    function [obj, opt, trainErr] = trainCmaes(obj, linear_hyp, X, y, lb, ub, f);
      % train with CMA-ES
      %
      global modelTrainNErrors;

      opt = []; fval = Inf; trainErr = false;
      cmaesopt.LBounds = lb';
      cmaesopt.UBounds = ub';
      if (length(obj.hyp.cov) > 2)
        % there is ARD covariance
        % try run cmaes for 500 funevals to get bounds for covariances
        MAX_DIFF = 2.5;
        cmaesopt.MaxFunEvals = 500;
        cmaesopt.SaveVariables = false;
        modelTrainNErrors = 0;
        try
          [opt, fval] = s_cmaes(f, linear_hyp', [0.3*(ub(1:(end-1)) - lb(1:(end-1))) 100]', cmaesopt);
        catch err
          fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
          trainErr = true;
          obj.nErrors = modelTrainNErrors;
          obj.trainGeneration = -1;
          return;
        end
        cov_median = median(opt(1:(end-4)));
        ub(1:(end-4)) = cov_median + MAX_DIFF;
        lb(1:(end-4)) = cov_median - MAX_DIFF;
        cmaesopt.LBounds = lb';
        cmaesopt.UBounds = ub';
      end
      cmaesopt.MaxFunEvals = 2000;
      try
        modelTrainNErrors = 0;
        [opt, fval] = s_cmaes(f, linear_hyp', [0.3*(ub(1:(end-1)) - lb(1:(end-1))) 100]', cmaesopt);
      catch err
        fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
        trainErr = true;
        obj.nErrors = modelTrainNErrors;
        obj.trainGeneration = -1;
        return;
      end
      obj.nErrors = modelTrainNErrors;
      obj.trainLikelihood = fval;
    end


    function [opts, nonlnc] = defaultFminconOpts(obj)
      % return the optimization parameters for fmincon()
      %
      opts = optimset('fmincon');
      opts = optimset(opts, ...
        'GradObj', 'on', ...
        'TolFun', 1e-6, ...
        'TolX', 1e-7, ...
        'MaxIter', 500, ...
        'MaxFunEvals', 500, ...
        'Display', 'off' ...
        );
      covarianceDim = length(obj.hyp.cov) - 1;
      if (covarianceDim > 1)
        % ARD
        opts = optimset(opts, 'Algorithm', 'interior-point');
        nonlnc = @nonlincons;
      else
        % ISOtropic
        opts = optimset(opts, 'Algorithm', 'trust-region-reflective');
        nonlnc = [];
      end
    end

    function [lb_hyp, ub_hyp] = defaultLUBounds(obj)
      % return default lower/upper bounds for GP model hyperparameter training
      %
      lb_hyp.cov = -2 * ones(size(obj.hyp.cov));
      lb_hyp.inf = log(1e-6);
      lb_hyp.lik = log(1e-6);
      lb_hyp.mean = -Inf;
      ub_hyp.cov = 25 * ones(size(obj.hyp.cov));
      ub_hyp.inf = log(10);
      ub_hyp.lik = log(10);
      ub_hyp.mean = Inf;
    end
  end
end

function [post, nlZ, dnlZ] = infExactCountErrors(hyp, mean, cov, lik, x, y)
  global modelTrainNErrors;
  try
    [post, nlZ, dnlZ] = infExact(hyp, mean, cov, lik, x, y);
  catch err
    modelTrainNErrors = modelTrainNErrors + 1;
    throw(err);
    % if (modelTrainNErrors > 20)
    %   throw(err);
    % end
  end
end

function [nlZ, dnlZ] = linear_gp(linear_hyp, s_hyp, inf, mean, cov, lik, x, y)
  hyp = rewrap(s_hyp, linear_hyp');
  [nlZ, s_dnlZ] = gp(hyp, inf, mean, cov, lik, x, y);
  dnlZ = unwrap(s_dnlZ)';
end

function [c, ceq] = nonlincons(x)
  % checks if the values x(1:(end-4)) are within 2.5 off median
  MAX_DIFF = 2.5;
  ceq = [];
  assert(size(x,2) == 1, 'Argument for nonlincons is not a vector');
  c = zeros(size(x));
  % test only for covariance parameters
  % TODO: are there always 4 more parameters?!
  c = abs(x(1:end-4) - median(x(1:end-4))) - MAX_DIFF;
end
