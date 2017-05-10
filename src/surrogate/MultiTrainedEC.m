classdef MultiTrainedEC < EvolutionControl
% The strategy of change is more strict than Kern's (2006) strategy:
% 
% if iters > 1
%   nInit = max(lambda, nInit + 2*n_b)
% elseif iters < 1
%   nInit = min(0.33, nInit - n*n_b)
% 
% where non-integer nInit is taken as probabilistic number between two
% neighbouring integers -- e.g. 2.33 is taken as 2 with probability 66%
% and 3 with probability 33%

  properties 
    model
    counteval

    nOrigInit
    rankFunc
    rankErrorThresh
    lastModel
    nTrainErrors
    maxTrainErrors
  end

  methods 
    function obj = MultiTrainedEC(surrogateOpts)
    % constructor
      obj.model = [];
      obj.counteval = 0;
      obj.nOrigInit = defopts(surrogateOpts, 'evoControlNOrigInit', 1);
      obj.rankFunc = defopts(surrogateOpts, 'evoControlRankFunc', @errRankMuOnly);
      obj.rankErrorThresh = defopts(surrogateOpts, 'evoControlRankErrorThresh', 0.1);
      obj.lastModel = [];
      obj.nTrainErrors = 0;
      obj.maxTrainErrors = defopts(surrogateOpts, 'evoControlMaxTrainErrors', 2);
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of multi-trained evolution control

      fitness_raw = [];
      arx = [];
      arxvalid = [];
      arz = [];
      surrogateStats = NaN(1, 8);
      nInit = obj.getProbNOrigInit();
      nOrigEvaled = 0;

      % extract cmaes state variables
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      BD = cmaesState.BD;
      dim = cmaesState.dim;
      mu = cmaesState.mu;
      countiter = cmaesState.countiter;

      yFinal = NaN(1,lambda);

      obj.model = ModelFactory.createModel(surrogateOpts.modelType, surrogateOpts.modelOpts, xmean');

      if (isempty(obj.model))
        % model could not be created :( use the standard CMA-ES
        return;
      end

      % sample lambda new points
      [xPop, xPopValid, zPop] = ...
          sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
      origEvaled = false(1,lambda);

      % get the training data from 'archive'
      nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);
      minTrainSize = obj.model.getNTrainData();
      [xTrain, yTrain] = archive.getClosestDataFromPoints(nArchivePoints, xPopValid', sigma, BD);
      if (size(yTrain, 1) < minTrainSize)
        % We don't have enough data for model training
        return;
      end

      % train the model 
      [obj, ok] = obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts);
      if (~ok)
        return;
      end

      [yModel1, sd2Model1] = obj.model.predict(xPopValid');

      % find the ordering of the points with highest expected ranking
      % error
      perm = expectedRankDiff(yModel1, sd2Model1, cmaesState.mu, obj.rankFunc);

      reevalID = false(1, lambda);
      reevalID(perm(1:nInit)) = true;
      xToReeval = xPop(:, reevalID);
      xToReevalValid = xPopValid(:, reevalID);
      zToReeval = zPop(:, reevalID);

      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, counteval] = ...
          sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, sigma, nInit, counteval, cmaesState, sampleOpts, varargin{:});
      fprintf('counteval: %d\n', counteval)
      yFinal(reevalID) = yNew;
      origEvaled(reevalID) = true; nOrigEvaled = sum(origEvaled);
      % update the Archive
      archive = archive.save(xNewValid', yNew', countiter);
      % the obj.models' dataset will be supplemented with this
      % new points during the next training using all the xTrain

      % retrain model
      % [xTrain, yTrain] = archive.getDataNearPoint(nArchivePoints, ...
      %     xmean', surrogateOpts.evoControlTrainRange, sigma, BD);
      [xTrain, yTrain] = archive.getClosestDataFromPoints(nArchivePoints, xPopValid(:,~origEvaled)', sigma, BD);

      [obj, ok] = obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts);
      if (~ok)
        return;
      end

      % predict with the retrained model and calculate 
      % change in ranking (estimate model error)
      [yModel2, sd2Model2] = obj.model.predict(xPopValid');
      [~, sort1] = sort(yModel1);
      ranking2   = ranking(yModel2);
      err = obj.rankFunc(ranking2(sort1), mu);
      lastErr = err;
      % Debug:
      fprintf('Ranking error: %f\n', err);

      iters = 0;
      n_b = 1;
      % while there is some non-trivial change in ranking, re-evaluate new poits
      while ((nOrigEvaled < lambda) && (err >= obj.rankErrorThresh))
        reevalID = false(1, lambda);
        iters = iters + 1;
        % find the ordering of the points with highest expected ranking error
        perm = expectedRankDiff(yModel2, sd2Model2, cmaesState.mu, obj.rankFunc);
        % do not re-evaluate what has already been evaluated
        perm(origEvaled) = [];
        % take the most interesting point
        reevalID(perm(1)) = true;
        xToReeval = xPop(:, reevalID);
        xToReevalValid = xPopValid(:, reevalID);
        zToReeval = zPop(:, reevalID);
        % original-evaluate the chosen one point
        [yNew, xNew, xNewValid, zNew, counteval] = ...
            sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, sigma, 1, counteval, cmaesState, sampleOpts, varargin{:});
        fprintf('counteval: %d\n', counteval)
        yFinal(reevalID) = yNew;
        origEvaled(reevalID) = true; nOrigEvaled = sum(origEvaled);
        % update the Archive
        archive = archive.save(xNewValid', yNew', countiter);

        % retrain model
        [xTrain, yTrain] = archive.getClosestDataFromPoints(nArchivePoints, xPopValid(:,~origEvaled)', sigma, BD);
        [obj, ok] = obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts);
        if (~ok)
          % training unsuccessful, raising nOrigInit
          obj.nOrigInit = min(lambda, obj.nOrigInit + n_b);
          return;
        end

        % predict with the retrained model and calculate 
        % change in ranking (estimate model error)
        [yModel3, sd2Model3] = obj.model.predict(xPopValid');
        [~, sort2] = sort(yModel2);
        ranking3   = ranking(yModel3);
        lastErr = err;
        err = obj.rankFunc(ranking3(sort2), mu);
        % Debug:
        fprintf('Ranking error (while cycle): %f\n', err);
        yModel2 = yModel3;
        sd2Model2 = sd2Model3;
      end

      if iters > 1
        % Debug:
        fprintf('nOrigInit = %.2f (was: %.2f)\n', min(lambda, obj.nOrigInit + 2*n_b), obj.nOrigInit);

        obj.nOrigInit = min(lambda, obj.nOrigInit + 2*n_b);
      elseif iters < 1
        % Debug:
        fprintf('nOrigInit = %.2f (was: %.2f)\n', max(0.33, obj.nOrigInit - n_b), obj.nOrigInit);

        obj.nOrigInit = max(0.33, obj.nOrigInit - n_b);
      end

      %{
      kendallErr = kendall(yPredict, yNew', 'type', 'Kendall');
      rmse = sqrt(sum((yPredict' - yNew).^2))/length(yNew);
      fprintf('  model-gener.: reevaluated %d pts, test RMSE = %f, Kendl. corr = %f.\n', obj.nOrigInit, rmse, kendallErr);
      surrogateStats = [rmse, kendallErr];
      %}

      if ~all(origEvaled)
        yModel = obj.model.predict((xPopValid(:, ~origEvaled))');
        % save also the last measured errRankMu
        surrogateStats = getModelStatistics(obj.model, cmaesState, surrogateOpts, sampleOpts, counteval, lastErr);

        yFinal(~origEvaled) = yModel';
        xNew = [xNew, xPop(:, ~origEvaled)];
        xNewValid = [xNewValid, xPopValid(:, ~origEvaled)];
        zNew = [zNew, zPop(:, ~origEvaled)];

        % shift the f-values:
        %   if the model predictions are better than the best original value
        %   in the model's dataset, shift ALL (!) function values
        %   Note: - all values have to be shifted in order to preserve predicted
        %           ordering of values
        %         - small constant is added because of the rounding errors
        %           when numbers of different orders of magnitude are summed
        fminDataset = min(obj.model.getDataset_y());
        fminModel = min(yModel);
        diff = max(fminDataset - fminModel, 0);
        yFinal = yFinal + 1.000001*diff;
      end

      % save the resulting re-evaluated population as the returning parameters
      fitness_raw = yFinal;
      arx = xPop;
      arxvalid = xPopValid;
      arz = zPop;

    end


    function nInit = getProbNOrigInit(obj)
      fracN = obj.nOrigInit - floor(obj.nOrigInit);
      plus = 0;
      if (fracN > 0)
        plus = (rand() < fracN);
      end
      nInit = floor(obj.nOrigInit) + plus;
    end

    function [obj, succ] = trainModelOrUseLast(obj, xTrain, yTrain, cmaesState, sampleOpts)
      succ = false;
      m = obj.model.train(xTrain, yTrain, cmaesState, sampleOpts);
      if (~m.isTrained())
        % Failure in training
        fprintf('Model cannot be trained. Trying last model...\n');
        if (~isempty(obj.lastModel) && obj.lastModel.isTrained() && obj.nTrainErrors < obj.maxTrainErrors)
          fprintf('Last model seems OK. Using it.\n');
          obj.model = obj.lastModel;
          obj.nTrainErrors = obj.nTrainErrors + 1;
          succ = true;
          return;
        else
          fprintf('Error. Last model cannot be used. Re-setting to normal CMA-ES.\n');
          return;
        end
      else
        % Successful training
        obj.nTrainErrors = 0;
        if (~isempty(obj.model) && obj.model.isTrained())
          obj.lastModel = obj.model;
        end
        obj.model = m;
        succ = true;
        return;
      end
    end
  end

end

function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
