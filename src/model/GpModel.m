classdef GpModel < Model & BayesianICModel
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
    shiftX
    options
    hyp
    nHyp
    covBounds
    likBounds
    meanFcn
    covFcn
    likFcn
    infFcn
    prior
    nErrors
    trainLikelihood
    cmaesCheckBounds
    nRestarts

    % mcmc-related options
    mcmcResults
    mcmcNSimu
    mcmcBurnin
    mcmcNChains
    nSimuPost
    predictFullPost

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

      % Training options
      %   Optimization Toolbox check
      obj.options.trainAlgorithm = defopts(obj.options, 'trainAlgorithm', 'fmincon');

      if (strcmpi(obj.options.trainAlgorithm, 'fmincon') ...
          && ~license('checkout', 'optimization_toolbox'))
        warning('GpModel: Optimization Toolbox license not available. Switching to minimize().');
        obj.options.trainAlgorithm = 'minimize';
      end

      obj.nRestarts = defopts(obj.options, 'nRestarts', 1);

      % GP hyper-parameter settings
      if (isfield(obj.options, 'hypOptions'))
        % if covariance function, starting points and bounds are all defined in a single struct, expand it
        if (isfield(obj.options.hypOptions, 'hyp'))
          obj.options.hyp = obj.options.hypOptions.hyp;
        end

        if (isfield(obj.options.hypOptions, 'covBounds'))
          obj.options.covBounds = obj.options.hypOptions.covBounds;
        end

        if (isfield(obj.options.hypOptions, 'covFcn'))
          obj.options.covFcn = obj.options.hypOptions.covFcn;
        end

        if (isfield(obj.options.hypOptions, 'likPrior'))
          obj.options.covPrior = obj.options.hypOptions.likPrior;
        end

        if (isfield(obj.options.hypOptions, 'covPrior'))
          obj.options.covPrior = obj.options.hypOptions.covPrior;
        end
      end

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
      if (iscell(obj.covFcn))
        covfcn = obj.covFcn{1};
      else
        covfcn = obj.covFcn;
      end

      prior = defopts(obj.options, 'prior', struct());
      fns = fieldnames(prior);
      for fi = 1:numel(fns)
        fn = fns{fi};
        if ischar(prior.(fn)) && exist(prior.(fn), 'var')
          prior.(fn) = str2func(prior.(fn));
        elseif iscell(prior.(fn))
          prior.(fn) = cellfun(@(x) myeval(x), prior.(fn), 'UniformOutput', false);
        else
          prior.(fn) = myeval(prior.(fn));
        end
      end
      obj.prior = prior;

      % expand covariance lengthscale hyperparameter according to
      % the dimension if ARD covariance specified and lengthscale is scalar
      if (length(obj.hyp.cov) >= 2 && (isequal(covfcn, @covSEard) ...
           || isequal(covfcn, @covMaternard) || isequal(covfcn, @covRQard) ...
           || isequal(covfcn, @covPERard) ...
         ))
        obj.hyp.cov = [repmat(obj.hyp.cov(1), obj.dim, 1); obj.hyp.cov(2:end)];
      end

      obj.meanFcn = str2func(defopts(obj.options, 'meanFcn', 'meanConst'));
      obj.likFcn  = str2func(defopts(obj.options, 'likFcn',  'likGauss'));
      obj.infFcn  = str2func(defopts(obj.options, 'infFcn',  'infExactCountErrors'));
      obj.options.normalizeY = defopts(obj.options, 'normalizeY', true);
      obj.options.centerX = defopts(obj.options, 'centerX', false);

      % GP hyperparameter bounds
      obj.covBounds = defopts(obj.options, 'covBounds', ...
          [-2*ones(size(obj.hyp.cov)), 25*ones(size(obj.hyp.cov))]);
      % expand also covariance Bounds if they do not respect ARD covariance
      if ((size(obj.covBounds,1) >= 2) && (isequal(covfcn, @covSEard) ...
          || isequal(covfcn, @covMaternard) || isequal(covfcn, @covRQard) ...
          || isequal(covfcn, @covPERard) ...
         ))
        obj.covBounds = [repmat(obj.covBounds(1,:), obj.dim, 1); obj.covBounds(2:end,:)];
      end

%       if (isequal(covfcn, @covADD))
%         dgs = obj.covFcn{2}{1};
%         dgs = dgs(dgs <= obj.dim);
%         obj.covFcn{2}{1} = dgs;
%         r = length(dgs);
% 
%         s = zeros(1, r);
%         for i = 1:r
%           % (d over dgs(i))
%           s(i) = prod(arrayfun(@(j) obj.dim - j, 0:(dgs(i) - 1))) / factorial(dgs(i));
%         end
% 
%         obj.hyp.cov = [repmat(obj.hyp.cov(1), obj.dim, 1); ...
%           log(obj.hyp.cov{2} ./ s')];
%         obj.covBounds = [repmat(obj.covBounds(1, :), obj.dim, 1); ...
%                          [log(obj.covBounds(2, 1) * ones(r, 1)) ...
%                           log(obj.covBounds(2, 2) ./ s')] ...
%         ];
%       end

      obj.likBounds = defopts(obj.options, 'likBounds', log([1e-3, 10]));
      obj.cmaesCheckBounds = defopts(obj.options, 'cmaesCheckBounds', true);
      obj.nHyp = numel(unwrap(obj.hyp));

      % MCMC-related settings
      obj.mcmcNChains = myeval(defopts(obj.options, 'mcmcNChains', 2));
      obj.mcmcNSimu = myeval(defopts(obj.options, 'mcmcNSimu', 500));
      obj.mcmcBurnin = myeval(defopts(obj.options, 'mcmcBurnin', 500));
      obj.nSimuPost = myeval(defopts(obj.options, 'nSimuPost', 20));
      obj.mcmcResults = struct('results', {{}}, 'chain', [], 's2chain', [], 'sschain', []);
      obj.predictFullPost = defopts(obj.options, 'predictFullPost', false);

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

    function k = getNParams(obj)
      k = numel(unwrap(obj.hyp));
    end

    function n = getNData(obj)
      % returns the real number of data used for training
      if ~isfield(obj.dataset, 'X') || isempty(obj.dataset.X)
        n = 0;
      else
        n = size(obj.dataset.X, 1);
      end
    end

    function lik = getNegLogML(obj, varargin)
      lik = obj.trainLikelihood;
    end

    function lik = getNegLogEst(obj, varargin)
      lik = obj.trainLikelihood;
    end

    function [lppd, ppd] = getLogPredDens(obj, varargin)
      % Computes a vector of log pointwise predictive densities averaged over
      % posterior hyperparams and a matrix of log pointwise predictive
      % densities:
      %   lppd  = log [(1/S) * \sum_s p(y_i | \theta_s)] -- N x 1 vector
      %   ppd   = log p(y_i | \theta_s)                  -- N x S matrix

      if ~(isequal(obj.likFcn, @likGauss) || isequal(obj.likFcn, 'likGauss'))
        error('Pointwise predictive density not implemented for non-Gaussian likelihoods.');
      end

      if nargin > 1
        % the number of posterior simulations
        nSimu = varargin{1};
      else
        nSimu = obj.nSimuPost;
      end

      assert(isfield(obj.mcmcResults, 'chain') && ~isempty(obj.mcmcResults.chain));
      chain = obj.mcmcResults.chain;
      idx = randsample(1:size(chain, 1), min(size(chain, 1), nSimu));

      X = obj.getDataset_X();
      y = (obj.getDataset_y() - obj.shiftY) / obj.stdY;
      N = size(X, 1);
      assert(all(size(y) == [N, 1]));

      % multiplicative factor in Gaussian density
      ppd_mult = zeros(N, nSimu);
      % the exponential term in Gaussian density
      ppd_exp = zeros(N, nSimu);
      % log predictive densities
      ppd = zeros(N, nSimu);

      for i = 1:nSimu
        hyp_s = rewrap(obj.hyp, chain(idx(i), :));
        [~, ~, fmu, fs2] = gp(hyp_s, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, ...
          X, y, X);
        assert(all(size(fmu) == [N, 1]));

        % the predictive density is a Gaussian with mean fmu and variance
        % (fs2 + sn2)^2
        sn2 = exp(hyp_s.lik);
        s = fs2 + sn2;
        ppd_mult(:, i) = 1 ./ (nSimu * sqrt(2*pi) * s);
        z = -0.5 * ((y - fmu).^2 ./ s.^2);
        ppd_exp(:, i) = z;
        ppd(:, i) = -0.5 * log(2*pi) - 2 * log(s) + z;
      end

      lppd = logsumexp(ppd_exp, ppd_mult, 2);
    end

    function chains = getChains(obj)
      chains = obj.mcmcResults.chain;
    end

    function sample = getNegLogLPost(obj)
      sample = reshape(obj.mcmcResults.sschain, numel(obj.mcmcResults.sschain), 1);
    end

    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      % TODO: *write this* properly according to dimension and
      %       covariance function set in options
      nData = 3 * obj.dim;
    end

    function obj = trainModel(obj, X, y, xMean, generation)
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

      if (obj.options.centerX)
        obj.shiftX = mean(obj.dataset.X, 1);
      else
        obj.shiftX = zeros(1, obj.dim);
      end

      % set the mean hyperparameter if is needed
      if (isequal(obj.meanFcn, @meanConst))
        obj.hyp.mean = median(yTrain);
      elseif (isequal(obj.meanFcn, @meanLinear))
        obj.hyp.mean = median(yTrain) / obj.dim * ones(obj.dim,1);
      end
      
      % eval hyperparameters
      if ischar(obj.hyp.lik)
        obj.hyp.lik = myeval(obj.hyp.lik);
      end

      if iscell(obj.hyp.cov)
        hyp_evaled = zeros(length(obj.hyp.cov), 1);
        i = 1;
        while i <= length(obj.hyp.cov)
          s = myeval(obj.hyp.cov{i});
          k = length(s);
          hyp_evaled(i:i+(k-1)) = s;
          i = i + k;
        end
        obj.hyp.cov = hyp_evaled;
      end

      % wrap the starting point for hyperparameters inside corresponding bounds
      obj.hyp.cov = min(obj.covBounds(:, 2), max(obj.covBounds(:,1), obj.hyp.cov));
      obj.hyp.lik = min(obj.likBounds(2), max(obj.likBounds(1), obj.hyp.lik));

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

        if (obj.nRestarts > 1)
          multi_start_points = lhsdesign(obj.nRestarts - 1, length(linear_hyp), 'smooth', 'on');
          multi_start_points = bsxfun(@times, multi_start_points, exp(ub)-exp(lb));
          multi_start_points = bsxfun(@plus, multi_start_points, exp(lb));
          linear_hyp = [linear_hyp; log(multi_start_points)];
        end
        fprintf('Linear hyp: \n');
        disp(linear_hyp);

        trainErrs = false(obj.nRestarts);
        for i = 1:obj.nRestarts
          fprintf('%d / %d optimization trial\n', i, obj.nRestarts);

          % gp() with linearized version of the hyper-parameters
          f = @(par) linear_gp(par, obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, linear_hyp_start, const_hyp_idx);

          if (strcmpi(alg, 'fmincon'))
            [obj, opt, lik, trainErr] = obj.trainFmincon(linear_hyp(i, :), obj.getDataset_X(), yTrain, lb, ub, f);

            if (trainErr)
              fprintf('Trying CMA-ES...\n');
              alg = 'cmaes';
            end
          end
          if (strcmpi(alg, 'cmaes'))
            [obj, opt, lik, trainErr] = obj.trainCmaes(linear_hyp(i, :), obj.getDataset_X(), yTrain, lb, ub, f);
          end
          trainErrs(i) = trainErr;

          % update the optimum
          if (~trainErr && lik < obj.trainLikelihood)
            linear_hyp_start(~const_hyp_idx) = opt;
            opt = linear_hyp_start;
            obj.trainGeneration = generation;
            obj.hyp = rewrap(obj.hyp, opt);
            obj.trainLikelihood = lik;
          elseif (trainErr)
            % DEBUG OUTPUT:
            fprintf('Optimization trial failed.\n');
          end
        end % multistart loop

        % DEBUG OUTPUT:
        fprintf('.. model-training likelihood = %f\n', obj.trainLikelihood);
        % disp(obj.hyp);

        if (all(trainErrs))
          % DEBUG OUTPUT:
          obj.trainGeneration = -1;
          fprintf(2, '.. model is not successfully trained, likelihood = %f\n', obj.trainLikelihood);
        end
      elseif strcmpi(alg, 'mcmc')
        linear_hyp_start = unwrap(obj.hyp);
        const_hyp_idx = false(1, length(linear_hyp_start));

        % double negative marginal log likelihood
        silent = true; % suppress numerical errors in likelihood
        likfun = @(par, data) 2 * linear_gp(par', obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, ...
          data.xdata, data.ydata, linear_hyp_start, const_hyp_idx, silent);
        
        % joint hyperprior
        priorfun = @(par, ~, ~) 2 * obj.hyperPrior(par, obj.hyp, obj.prior);

        modelTrainNErrors = 0;

        for i = 1:obj.mcmcNChains
          if i > 1
            % add Gaussian perturbation
            linear_hyp_start1 = mvnrnd(linear_hyp_start, 0.5 * eye(obj.nHyp))';
          else
            linear_hyp_start1 = linear_hyp_start;
          end
          [results, chain, s2chain, sschain] = obj.trainMcmc(linear_hyp_start1, obj.getDataset_X(), yTrain, likfun, priorfun);
          obj.mcmcResults.results = [obj.mcmcResults.results results];
          obj.mcmcResults.chain = [obj.mcmcResults.chain; chain];
          obj.mcmcResults.s2chain = [obj.mcmcResults.s2chain; s2chain];
          obj.mcmcResults.sschain = [obj.mcmcResults.sschain; sschain];
        end

        obj.nErrors = modelTrainNErrors;

        % compute a Bayes estimate of hyperparameters from the chain
        % using an estimator function
        estfun = @median;
        expchain = exp(chain);
        expchain = expchain(~any(isinf(expchain), 2), :);
        est = log(feval(estfun, expchain));

        hyp_est = struct( ...
          'val', est, ...
          'lik', 0.5 * likfun(est, struct('xdata', obj.getDataset_X(), 'ydata', yTrain)), ...
          'prior', 0.5 * priorfun(est) ...
        );

        % a point estimate
        obj.hyp = rewrap(obj.hyp, hyp_est.val);
        obj.trainLikelihood = hyp_est.lik;

        if ~isinf(hyp_est.val)
          obj.trainGeneration = generation;
        else
          obj.trainGeneration = -1;
        end
      else
        error('GpModel.train(): train algorithm "%s" is not known.\n', alg);
      end
    end

    function [y, sd2] = modelPredict(obj, X)
      % predicts the function values in new points X
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        XWithShift = XWithShift - obj.shiftX; % centering Xs
        % prepare the training set (if was normalized for training)
        yTrain = (obj.getDataset_y() - obj.shiftY) / obj.stdY;

        if obj.predictFullPost
          C = size(obj.mcmcResults.chain, 1);
          idx = randsample(1:C, obj.nSimuPost);
          hypSimu = obj.mcmcResults.chain(idx, :);
          N = size(X, 1);
          mu = zeros(N, obj.nSimuPost);
          s2 = zeros(N, obj.nSimuPost);
          gpfail = true(1, obj.nSimuPost);
          for s = 1:obj.nSimuPost
            try
              h = rewrap(obj.hyp, hypSimu(s, :));
              [~, ~, fmu, fs2] = gp(h, obj.infFcn, obj.meanFcn, ...
                obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, XWithShift);

              % un-normalize in the f-space (if there is any)
              mu(:, s) = fmu .* obj.stdY + obj.shiftY;
              s2(:, s) = fs2 .* (obj.stdY)^2;
              gpfail(s) = false;
            catch
              gpfail(s) = true;
            end
          end

          % average the predictive mvns over the hyperposterior sample
          y = mean(mu(:, ~gpfail), 2);
          sd2 = mean(s2(:, ~gpfail), 2) / obj.nSimuPost;
        else
          % calculate GP models' prediction in X
          [y, gp_sd2] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, XWithShift);
          % un-normalize in the f-space (if there is any)
          y = y * obj.stdY + obj.shiftY;
          sd2 = gp_sd2 * (obj.stdY)^2;
        end

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
        fprintf(2, 'GpModel.predict(): the model is not yet trained!\n');
      end
    end

    function [y, sd2] = modelPredictFullPost(obj, X)
      % predicts the function values in new points X averaging over full posterior
      % p(f, \theta)
      if (obj.isTrained())
        % apply the shift if the model is already shifted
        XWithShift = X - repmat(obj.shiftMean, size(X,1), 1);
        XWithShift = XWithShift - obj.shiftX; % centering Xs
        % prepare the training set (if was normalized for training)
        yTrain = (obj.getDataset_y() - obj.shiftY) / obj.stdY;

        C = size(obj.mcmcResults.chain, 1);
        idx = randsample(1:C, obj.nSimuPost);
        hypSimu = obj.mcmcResults.chain(idx, :);
        hypBak = obj.hyp;

        N = size(X, 1);
        mu = zeros(obj.nSimuPost, N);
        s2 = zeros(obj.nSimuPost, obj.nSimuPost, N);
        for s = 1:obj.nSimuPost
          obj.hyp = rewrap(obj.hyp, hypSimu(s, :));
          [~, ~, fmu, fs2] = gp(obj.hyp, obj.infFcn, obj.meanFcn, obj.covFcn, obj.likFcn, obj.getDataset_X(), yTrain, XWithShift);

          % un-normalize in the f-space (if there is any)
          mu(s, :) = fmu * obj.stdY + obj.shiftY;
          s2(s, :) = fs2 * (obj.stdY)^2;
        end

        y = mean(mu, 2);
        sd2 = sum(mu, 2) / (obj.nSimuPost^2);

        obj.hyp = hypBak;
      else
        y = []; sd2 = [];
        fprintf(2, 'GpModel.predict(): the model is not yet trained!\n');
      end
    end

    function X = getDataset_X(obj)
      X = obj.dataset.X - obj.shiftX;
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

    function [obj, opt, fval, trainErr] = trainFmincon(obj, linear_hyp, X, y, lb, ub, f)
      % train with Matlab's fmincon() from the Optimization toolbox
      %
      global modelTrainNErrors;
      trainErr = false;
      opt = [];
      fval = [];

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
        fprintf('Model training (fmincon), init fval = %.2f\n', initial);
        try
          modelTrainNErrors = 0;
          [opt, fval] = fmincon(f, linear_hyp', [], [], [], [], lb, ub, [], fminconOpts);
          obj.nErrors = modelTrainNErrors;
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

    function [obj, opt, fval, trainErr] = trainCmaes(obj, linear_hyp, X, y, lb, ub, f)
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
      if (length(obj.hyp.cov) > 2 && obj.cmaesCheckBounds)
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
        return;
      end
      if (isnan(fval))
        % final likelihood is not a valid value, the model is probably
        % not trained, do not use it
        trainErr = true;
      end
      obj.nErrors = modelTrainNErrors;
    end


    function [results, chain, s2chain, sschain] = trainMcmc(obj, linear_hyp_start, X, y, likfun, priorfun)
      data.xdata = X;
      data.ydata = y;

      model = struct( ...
        'ssfun', likfun, ...
        'priorfun', priorfun ...
      );

      % parameter structure
      h = length(linear_hyp_start);
      params = cell(1, h);
      names = cell(1, h);

      j = 1;
      for name_cell = fieldnames(orderfields(obj.hyp))'
        name = name_cell{:};
        k = numel(obj.hyp.(name));
        names(j:j+k-1) = arrayfun(@(d) sprintf([name '%d'], d), 1:k, 'UniformOutput', false);
        j = j + k;
      end

      for i = 1:h
        params{i} = {names{i}, linear_hyp_start(i)};
      end

      mcmcOpts = struct( ...
        'nsimu',         obj.mcmcNSimu, ...
        'burnintime',    obj.mcmcBurnin, ...
        'verbosity',     0, ...
        'updatesigma',   true, ...
        'waitbar',       false ...
      );
      [results, chain, s2chain, sschain] = mcmcrun(model, data, params, mcmcOpts);
    end


    function [opts, nonlnc] = defaultFminconOpts(obj, lb, ub)
      % return the optimization parameters for fmincon()
      %
      opts = optimset('fmincon');
      opts = optimset(opts, ...
        'GradObj', 'on', ...
        'TolFun', 1e-6, ...
        'TolX', 1e-7, ...
        'MaxIter', 1000, ...
        'MaxFunEvals', 3000, ...
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

  methods (Static)
    function prob = logHyperPrior(linear_hyp, hyp)
      % A log product of hyperparameter priors.
      % All hyp.likPrior and hyp.covPrior{:} functions are assumed to be
      % densities.
      hyp_s = rewrap(hyp, exp(linear_hyp'));
      logp = zeros(1, 1 + numel(hyp.cov));

      assert(numel(hyp.lik) == 1);

      for i = 1:numel(hyp.covPrior)
        priorFcn = hyp.covPrior{i};
        logp(i) = log(priorFcn(hyp_s.cov(i)));
      end
      logp(end) = hyp.likPrior(hyp_s.lik);

      prob = sum(logp);
    end

    function nlp = hyperPrior(linear_hyp, hyp, prior)
      inf = @infZeros;
      hyp_s = vec2any(hyp, linear_hyp');

      % infPrior computes nlZ and distracts log prior value
      % using infZeros efficiently computes log prior, since supplied nlZ is zero
      % this hack preserves flexibility of gpml prior specifications,
      % but still evaluates prior derivatives, whis is not needed here
      % FIXME: write a proper evalPrior function working with gpml-style prior specs
      nlp = infPriorOnly(inf, prior, hyp_s);
    end

  end

end

function [post, nlZ, dnlZ] = infZeros(~, ~, ~, ~, ~, ~)
  post = []; nlZ = 0; dnlZ = 0;
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

function [nlZ, dnlZ] = linear_gp(linear_hyp, s_hyp, infFcn, mean, cov, lik, x, y, linear_hyp_start, const_hyp_idx, varargin)
  % extend the vector of parameters by constant elements
  % identified in the vector linear_hyp_start by indices const_hyp_idx
  linear_hyp_start(~const_hyp_idx) = linear_hyp;
  linear_hyp = linear_hyp_start;

  hyp = rewrap(s_hyp, linear_hyp');
 
  try
    if nargout <= 1
      % compute negative marginal log likelihood
      nlZ = gp(hyp, infFcn, mean, cov, lik, x, y);
    else
      % compute also negative marginal log likelihood derivatives
      [nlZ, s_dnlZ] = gp(hyp, infFcn, mean, cov, lik, x, y);
      dnlZ = unwrap(s_dnlZ)';
      dnlZ = dnlZ(~const_hyp_idx);
    end
  catch err
    if nargin > 10
      % fail silently (return inf)
      silentflag = varargin{1};
    else
      silentflag = false;
    end

    if silentflag
      if nargout <= 1
        nlZ = inf;
      else
        nlZ = inf;
        dnlZ = zeros(length(linear_hyp), 1);
        dnlZ = dnlZ(~const_hyp_idx);
      end
    else
      throw(err);
    end
  end
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
