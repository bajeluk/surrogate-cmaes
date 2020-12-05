classdef GpModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling

    % GpModel specific fields
    stdY                  % standard deviation of Y in training set, for normalizing output
    options
    hyp
    covBounds
    likBounds
    meanFcn
    covFcn
    likFcn
    infFcn
    nErrors
    trainLikelihood

    % Dimensionality-reduction specific fields
    dimReduction          % Reduce dimensionality for model by eigenvectors
                          % of covatiance matrix in percentage
    reductionMatrix       % Matrix used for dimensionality reduction
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
        fprintf('GpModel: Optimization Toolbox license not available. Switching to minimize().\n');
        obj.options.trainAlgorithm = 'minimize';
      end

      % GP hyper-parameter settings
      if (~isfield(obj.options, 'hyp'))
        obj.options.hyp = struct();
      end
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
      % expand covariance lengthscale hyperparameter according to
      % the dimension if ARD covariance specified and lengthscale is scalar
      if (iscell(obj.covFcn)) covfcn = obj.covFcn{1};
      else covfcn = obj.covFcn; end
      if (length(obj.hyp.cov == 2) && (isequal(covfcn, @covSEard) ...
          || isequal(covfcn, @covMaternard)))
        obj.hyp.cov = [obj.hyp.cov(1)*ones(obj.dim, 1); obj.hyp.cov(2)];
      end
      obj.meanFcn = str2func(defopts(obj.options, 'meanFcn', 'meanConst'));
      obj.likFcn  = str2func(defopts(obj.options, 'likFcn',  'likGauss'));
      obj.infFcn  = str2func(defopts(obj.options, 'infFcn',  'infExactCountErrors'));
      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);

      % GP hyperparameter bounds
      obj.covBounds = defopts(obj.options, 'covBounds', ...
          [-2*ones(size(obj.hyp.cov)), 25*ones(size(obj.hyp.cov))]);
      % expand also covariance Bounds if they do not respect ARD covariance
      if ((size(obj.covBounds,1) == 2) && (isequal(covfcn, @covSEard) ...
          || isequal(covfcn, @covMaternard)))
        obj.covBounds = [repmat(obj.covBounds(1,:), obj.dim, 1); obj.covBounds(2:end,:)];
      end
      obj.likBounds = defopts(obj.options, 'likBounds', log([1e-6, 10]));

      % wrap the starting point for hyperparameters inside corresponding bounds
      obj.hyp.cov = min(obj.covBounds(:, 2), max(obj.covBounds(:,1), obj.hyp.cov));
      obj.hyp.lik = min(obj.likBounds(2), max(obj.likBounds(1), obj.hyp.lik));

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.dimReduction = defopts(modelOptions, 'dimReduction', 1);      % 1.0 == no dimensionality reduction
    end

    function obj = clone(obj, obj2)
    % Take all fields except function handles from obj2
      fnames = fieldnames(obj2);
      for i = 1:length(fnames)
        ff = fnames{i};
        if (isempty(strfind(ff, 'Fcn')))
          obj.(ff) = obj2.(ff);
        end
      end
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = 3 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation, ~)
      % train the GP model based on the data (X,y)
      % TODO
      %   [ ] implement choosing the best covariance function according to
      %       the test ordinal regression capabilities
      global modelTrainNErrors;

      assert(size(xMean,1) == 1, '  GpModel.train(): xMean is not a row-vector.');
      obj.trainMean = xMean;
      if (~isempty(X) && ~isempty(y))
        obj.dataset.X = X;
        obj.dataset.y = y;
      end

      % normalize y if specified, @meanZero, or if large y-scale
      % (at least for CMA-ES hyperparameter optimization)
      if (~obj.options.normalizeY ...
          && (isequal(obj.meanFcn, @meanZero) || (max(obj.dataset.y) - min(obj.dataset.y)) > 1e4))
        fprintf(2, 'Y-Normalization is switched ON for @meanZero covariance function of large Y-scale.\n');
        obj.options.normalizeY = true;
      end
      if (obj.options.normalizeY)
        obj.shiftY = mean(obj.dataset.y);
        obj.stdY  = std(obj.dataset.y);
        yTrain = (obj.dataset.y - obj.shiftY) / obj.stdY;
      else
        obj.shiftY = 0;
        obj.stdY  = 1;
        yTrain = obj.dataset.y;
      end

      % set the mean hyperparameter if is needed
      if (isequal(obj.meanFcn, @meanConst))
        obj.hyp.mean = median(yTrain);
      elseif (isequal(obj.meanFcn, @meanLinear))
        obj.hyp.mean = median(yTrain) / obj.dim * ones(obj.dim,1);
      end
      

      alg = obj.options.trainAlgorithm;

      if (strcmpi(alg, 'minimize'))
        [obj, fval] = obj.trainMinimize(obj.getDataset_X(), yTrain);
        if (fval < Inf)
          obj.trainGeneration = generation;
        else
          obj.trainGeneration = -1;
        end

      elseif (strcmpi(alg, 'fmincon') ...
              || strcmp(alg, 'cmaes'))
        % lower and upper bounds
        [lb_hyp, ub_hyp] = obj.getLUBounds(yTrain, obj.hyp);
        lb = unwrap(lb_hyp)';
        ub = unwrap(ub_hyp)';
        opt = [];

        linear_hyp = unwrap(obj.hyp)';
        l_cov = length(obj.hyp.cov);

        % if some parameters are held constant
        const_hyp_idx = (lb == ub);
        linear_hyp_start = linear_hyp;
        linear_hyp = linear_hyp(~const_hyp_idx);
        lb = lb(~const_hyp_idx);
        ub = ub(~const_hyp_idx);

        % gp() with linearized version of the hyper-parameters
        f = @(par) linear_gp(par, obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, linear_hyp_start, const_hyp_idx);

        if (strcmpi(alg, 'fmincon'))
          [obj, opt, trainErr] = obj.trainFmincon(linear_hyp, obj.getDataset_X(), yTrain, lb, ub, f);

          if (trainErr)
            % fprintf('Trying CMA-ES...\n');
            alg = 'cmaes';
          end
        end
        if (strcmpi(alg, 'cmaes'))
          [obj, opt, trainErr] = obj.trainCmaes(linear_hyp, obj.getDataset_X(), yTrain, lb, ub, f);
          if (trainErr)
            % DEBUG OUTPUT:
            fprintf(2, '.. model is not successfully trained, likelihood = %f\n', obj.trainLikelihood);
            return;
          end
        end

        linear_hyp_start(~const_hyp_idx) = opt;
        opt = linear_hyp_start;
        obj.trainGeneration = generation;
        obj.hyp = rewrap(obj.hyp, opt);

        % DEBUG OUTPUT:
        % fprintf('.. model-training likelihood = %f\n', obj.trainLikelihood);
        % disp(obj.hyp);
      else
        error('GpModel.train(): train algorithm "%s" is not known.\n', alg);
      end
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        % prepare the training set (if was normalized for training)
        yTrain = (obj.getDataset_y() - obj.shiftY) / obj.stdY;
        % calculate GP models' prediction in X
        [y, gp_sd2] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, XWithShift);
        % un-normalize in the f-space (if there is any)
        y = y * obj.stdY + obj.shiftY;
        sd2 = gp_sd2 * (obj.stdY)^2;

        % % Calculate POI if it should be used
        % if (obj.options.usePOI)
        %   % return -POI , because the smaller y-value (bigger probability) the better
        %   y = - getPOI(X, y, dev, min(obj.getDataset_y()));
        %   dev = zeros(size(dev));
        % end
        %
        % % Calculate EI if it should be used
        % if (obj.options.useEI)
        %   % EI should be negative in promising regions, the lower the better
        %   y = getEI(X, y, dev, min(obj.getDataset_y()));
        %   dev = zeros(size(dev));
        % end
      else
        y = []; sd2 = [];
        fprintf(2, 'GpModel.(): the model is not yet trained!\n');
      end
    end
    
    function [x] = minimumX(obj, archive)
        ub = obj.sampleOpts.ubounds;
        lb = obj.sampleOpts.lbounds;
        
        cmaesopt.LBounds = lb;
        cmaesopt.UBounds = ub;
        cmaesopt.SaveVariables = false;
        cmaesopt.LogModulo = 0;
        cmaesopt.DispModulo = 0;
        cmaesopt.DispFinal = 0;
        cmaesopt.Seed = 'inherit';
        sigma = [0.3*(ub - lb)];
        % sigma(end) = min(10*mean(sigma(1:end-1)), sigma(end));
        % there is ARD covariance
        % try run cmaes for 500 funevals to get bounds for covariances
        MAX_DIFF = 2.5;
        cmaesopt.MaxFunEvals = 500;
        modelTrainNErrors = 0;
        
        
        eval_func = @(X) obj.predict(X');
        [opt, fval] = s_cmaes(eval_func, obj.trainMean, sigma, cmaesopt);        
        x = opt';
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
      % fprintf('  ... minimize() %f --> %f in %d iterations.\n', fval(1), fval(end), iters);
      warning('on');

      obj.nErrors = modelTrainNErrors;
      obj.trainLikelihood = fval(end);
      obj.hyp = hyp_;
    end

    function [obj, opt, trainErr] = trainFmincon(obj, linear_hyp, X, y, lb, ub, f)
      % train with Matlab's fmincon() from the Optimization toolbox
      %
      global modelTrainNErrors;
      trainErr = false;
      opt = [];

      [fminconOpts, nonlnc] = obj.defaultFminconOpts(lb, ub);
      try
        initial = f(linear_hyp');
      catch err
        initial = NaN;
      end
      if isnan(initial)
        % the initial point is not valid
        fprintf('  GpModel.train(): fmincon -- initial point is not valid.\n');
        trainErr = true;
      else
        % training itself
        % fprintf('Model training (fmincon), init fval = %f\n', num2str(initial));
        try
          modelTrainNErrors = 0;
          [opt, fval] = fmincon(f, linear_hyp', [], [], [], [], lb, ub, nonlnc, fminconOpts);
          obj.nErrors = modelTrainNErrors;
          obj.trainLikelihood = fval;
          if (isnan(fval)  ||  initial - fval < 0.1)
            % final likelihood is not a valid value or
            % the shift in likelihood is almost none, the model is probably
            % not trained, do not use it
            trainErr = true;
          end
        catch err
          obj.nErrors = modelTrainNErrors;
          fprintf(2, '  GpModel.train() ERROR: fmincon() ended with an exception: %s\n', err.message);
          trainErr = true;
        end
      end
    end

    function [obj, opt, trainErr] = trainCmaes(obj, linear_hyp, X, y, lb, ub, f)
      % train with CMA-ES
      %
      global modelTrainNErrors;

      opt = []; fval = Inf; trainErr = false;
      cmaesopt.LBounds = lb';
      cmaesopt.UBounds = ub';
      cmaesopt.SaveVariables = false;
      cmaesopt.LogModulo = 0;
      cmaesopt.DispModulo = 0;
      cmaesopt.DispFinal = 0;
      cmaesopt.Seed = 'inherit';
      sigma = [0.3*(ub - lb)]';
      % sigma(end) = min(10*mean(sigma(1:end-1)), sigma(end));
      if (length(obj.hyp.cov) > 2)
        % there is ARD covariance
        % try run cmaes for 500 funevals to get bounds for covariances
        MAX_DIFF = 2.5;
        cmaesopt.MaxFunEvals = 500;
        modelTrainNErrors = 0;
        try
          [opt, fval] = s_cmaes(f, linear_hyp', sigma, cmaesopt);
        catch err
          fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
          trainErr = true;
          obj.nErrors = modelTrainNErrors;
          obj.trainGeneration = -1;
          return;
        end
        cov_median = median(opt(1:obj.dim));
        % ub(1:obj.dim) = cov_median + MAX_DIFF;
        ub(1:obj.dim) = min(max(opt(1:obj.dim)', linear_hyp(1:obj.dim)) + MAX_DIFF, ub(1:obj.dim));
        % lb(1:obj.dim) = cov_median - MAX_DIFF;
        lb(1:obj.dim) = max(min(opt(1:obj.dim)', linear_hyp(1:obj.dim)) - MAX_DIFF, lb(1:obj.dim));
        cmaesopt.LBounds = lb';
        cmaesopt.UBounds = ub';
        sigma(1:obj.dim) = [0.3*(ub(1:obj.dim) - lb(1:obj.dim))]';
      end
      cmaesopt.MaxFunEvals = 2000;
      try
        modelTrainNErrors = 0;
        [opt, fval] = s_cmaes(f, linear_hyp', sigma, cmaesopt);
      catch err
        fprintf(2, 'GpModel.train() ERROR: CMA-ES ended with an exception: %s\n', err.message);
        trainErr = true;
        obj.nErrors = modelTrainNErrors;
        obj.trainGeneration = -1;
        return;
      end
      if (isnan(fval))
        % final likelihood is not a valid value, the model is probably
        % not trained, do not use it
        trainErr = true;
      end
      obj.nErrors = modelTrainNErrors;
      obj.trainLikelihood = fval;
    end


    function [opts, nonlnc] = defaultFminconOpts(obj, lb, ub)
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
        % ARD or a parameter with fixed value
        opts = optimset(opts, 'Algorithm', 'interior-point');
        nonlnc = @nonlincons;
      else
        % ISOtropic
        opts = optimset(opts, 'Algorithm', 'trust-region-reflective');
        nonlnc = [];
      end
    end

    function [lb_hyp, ub_hyp] = getLUBounds(obj, yTrain, startHyp)
      % return lower/upper bounds for GP model hyperparameter training
      %
      lb_hyp.cov = obj.covBounds(:,1);
      ub_hyp.cov = obj.covBounds(:,2);
      lb_hyp.lik = obj.likBounds(1);
      ub_hyp.lik = obj.likBounds(2);
      % set bounds for mean hyperparameter
      if (isequal(obj.meanFcn, @meanConst))
        minY = min(yTrain);
        maxY = max(yTrain);
        lb_hyp.mean = minY - 2*(maxY - minY);
        ub_hyp.mean = minY + 2*(maxY - minY);
      elseif (isequal(obj.meanFcn, @meanLinear))
        min_y = min(yTrain);
        max_y = max(yTrain);
        lb_hyp.mean = zeros(size(startHyp.mean));
        ub_hyp.mean = zeros(size(startHyp.mean));
        for i=1:obj.dim
          % max_x -- max of each dimension from dataset_X
          dataset_X = obj.getDataset_X();
          max_x = max(dataset_X(:,i));
          min_x = min(dataset_X(:,i));
          max_tg = (max_y - min_y) / (max_x - min_x);
          lb_hyp.mean(i) = -5 * max_tg;
          ub_hyp.mean(i) = 5 * max_tg;
        end
        lb_hyp.mean = min(lb_hyp.mean, startHyp.mean);
        ub_hyp.mean = max(ub_hyp.mean, startHyp.mean);
      end
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

function [nlZ, dnlZ] = linear_gp(linear_hyp, s_hyp, inf, mean, cov, lik, x, y, linear_hyp_start, const_hyp_idx)
  % extend the vector of parameters by constant (i.e. not optimized) elements
  % taken from the vector of initial values
  linear_hyp_start(~const_hyp_idx) = linear_hyp;
  linear_hyp = linear_hyp_start;

  hyp = rewrap(s_hyp, linear_hyp');
  [nlZ, s_dnlZ] = gp(hyp, inf, mean, cov, lik, x, y);
  dnlZ = unwrap(s_dnlZ)';
  dnlZ = dnlZ(~const_hyp_idx);
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
