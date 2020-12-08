classdef ModelAssistedEC < IndividualEC
  % Implements Model Assisted ES
  %
  % Reference: H. Ulmer, F. Streichert, A. Zell,
  % "Evolution Strategies assisted by Gaussian Processes with improved
  % Pre-Selection Criterion",
  % IEEE Congress on Evolutionary Computation,CEC 2003: 692-699

  methods
    function obj = ModelAssistedEC()
      % constructor
      obj@IndividualEC();
    end

    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
      % Run one generation of Model Assisted ES's evolution control

      fitness_raw = [];
      arx = [];
      arxvalid = [];
      arz = [];
      surrogateStats = NaN(1, 2);

      % extract cmaes state variables
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      dim = cmaesState.dim;
      countiter = cmaesState.countiter;

      % create model
      obj.model = ModelFactory.createModel(surrogateOpts.modelType, surrogateOpts.modelOpts, xmean');

      if (isempty(obj.model))
        % model could not be created :( use the standard CMA-ES
        warning('ModelAssistedEC.runGeneration(): could not create a model');
        return;
      end

      if isfield(surrogateOpts.modelOpts, 'trainsetSizeMax')
        nArchivePoints = myeval(surrogateOpts.modelOpts.trainsetSizeMax);
      else
        nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);
      end

      [xTrain, yTrain] = ModelAssistedEC.getRecentData(archive, nArchivePoints);

      if (size(xTrain, 1) < nArchivePoints)
        warning('not enough points (%d) for training (required %d)', size(xTrain, 1), nArchivePoints);
        return;
      end

      obj.model = obj.model.train(xTrain, yTrain, cmaesState, sampleOpts);

      if (~obj.model.isTrained())
        warning('ModelAssistedEC.runGeneration(): model not trained');
        return;
      end

      % sample the enlarged population of size 'gamma * lambda'
      extendSize = ceil(surrogateOpts.evoControlIndividualExtension ...
          * lambda);
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(sigma, extendSize, cmaesState, sampleOpts);

      % calculate the model prediction for the extended population
      yExtend = obj.getModelOutput(xExtend');

      % choose lambda best points
      [xToReeval, xToReevalValid, zToReeval] = ...
          SurrogateSelector.choosePointsToReevaluate(...
          xExtend, xExtendValid, zExtend, yExtend, lambda, 0);

      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, counteval] = ...
          sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, sigma, lambda, counteval, cmaesState, sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
      fprintf('counteval: %d\n', counteval)
      % update the Archive
      archive = archive.save(xNewValid', yNew', countiter);
      % the obj.models' dataset will be supplemented with this
      % new points during the next training using all the xTrain
      % calculate the models' precision
      yPredict = obj.model.predict(xNewValid');
      kendall = corr(yPredict, yNew', 'type', 'Kendall');
      rmse = sqrt(sum((yPredict' - yNew).^2))/length(yNew);
      fprintf('  model-gener.: %d preSamples, reevaluated %d pts, test RMSE = %f, Kendl. corr = %f.\n', 0, lambda, rmse, kendall);
      surrogateStats = [rmse kendall];

      % save the resulting re-evaluated population as the returning parameters
      fitness_raw = [fitness_raw yNew];
      arx = [arx xNew];
      arxvalid = [arxvalid xNewValid];
      arz = [arz zNew];
    end % function

    function [output, y] = getModelOutput(obj, X)
      % overrides 'Model.getModelOutput' especially in the 'poi' branch
      % where 'target' is set exactly to 'fmin'
      [y, sd2] = obj.model.predict(X);
      fmin = min(obj.model.dataset.y);

      switch lower(obj.model.predictionType)
        case 'fvalues' % mean function values
          output = y;

        case 'poi' % probability of improvement
          target = fmin;
          output = getPOI(X, y, sd2, target);
      end
    end % function
  end % methods

  methods (Static)
    function [X, y] = getRecentData(archive, n)
      % Get up to 'n' most recently added points from 'archive'.
      nData = length(archive.y);
      X = []; y = [];

      if (nData == 0)
        return;
      else
        X = archive.X(max(1, nData - n + 1):nData, :);
        y = archive.y(max(1, nData - n + 1):nData, :);
        return;
      end
    end % function

    function Kmn = covFcn(Xm, Xn, theta)
      % A covariance function defined as sum of SEard, constant
      % shift and a white noise:
      %
      % Kmn(x_i, x_j) = sf2 * exp(-(x_i - x_q)'*diag(1./ell)*(x_i - x_j)/2)
      %                 + sc2
      %                 + sn2*delta(x_i, x_j)
      %
      % where ell = (ell_1^2,...ell_D^2) is vector of ARD parameters, D is
      % dimension of input space, sf2 is signal variance, sc2 is constant shift
      % and sn2 is added white noise.
      %
      % Theta is a vector of positive numbers:
      % theta = [ log(ell_1)
      %           log(ell_2)
      %            .
      %           log(ell_D)
      %           log(sqrt(sf2))
      %           log(sqrt(sc2))
      %           log(sqrt(sn2)) ]
      %
      assert(size(Xm, 2) == size(Xn, 2), 'Dimensions must agree');
      dim = size(Xm, 2);

      % tolerance on squared distance of vectors when their considered equal
      tol = 1e-9;

      assert(all(size(theta) == [dim+3 1]), ['theta must be a row vector of ''' dim + 3 ''' hyperparameters']);

      ell = exp(theta(1:dim));
      sf2 = exp(2*theta(dim+1));
      sc2 = exp(2*theta(dim+2));
      sn2 = exp(2*theta(dim+3));

      Kmn = sq_dist(diag(1./ell) * Xm', diag(1./ell) * Xn');
      Kmn = sf2 * exp(-Kmn/2);
      Kmn = bsxfun(@plus, Kmn, sc2);

      delta = bsxfun(@lt, sq_dist(Xm', Xn'), tol^2);
      Kmn(delta) = Kmn(delta) + sn2;
    end
  end % methods (Static)
end % classdef

function res = myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
