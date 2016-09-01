classdef MultiTrainedEC < EvolutionControl
  properties 
    model
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
      obj.nOrigInit = defopts(surrogateOpts, 'evoControlNOrigInit', 1);
      obj.rankFunc = defopts(surrogateOpts, 'evoControlRankFunc', @errRankMuOnly);
      obj.rankErrorThresh = defopts(surrogateOpts, 'evoControlRankErrorThresh', 0.1);
      obj.lastModel = [];
      obj.nTrainErrors = 0;
      obj.maxTrainErrors = defopts(surrogateOpts, 'evoControlMaxTrainErrors', 2);
    end

    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
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
      if (~obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts))
        return;
      end

      [yModel1, sd2Model1] = obj.model.predict(xPopValid');

      % find the ordering of the points with highest expected ranking
      % error
      perm = obj.mostProbableRankDiff(yModel1, sd2Model1, cmaesState.mu);

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

      if (~obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts))
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
        perm = obj.mostProbableRankDiff(yModel2, sd2Model2, cmaesState.mu);
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
        if (~obj.trainModelOrUseLast(xTrain, yTrain, cmaesState, sampleOpts))
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
        fminDataset = min(obj.model.dataset.y);
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

    function perm = mostProbableRankDiff(obj, yPredict, sd2Predict, mu)
      %MOSTPROBABLERANKDIFF returns permutation of points which causes highest expected difference in ranking
      %
      % y       mean predicions returned by the GP/rforest model
      % sd2     predicited variances for respective predictions in y
      % mu      CMA-ES' mu parameter -- how many points are taken for
      %         updating mean and rank-mu update of covariance matrix
      %
      % perm    Returns permutation of the points in the order from maximum
      %         expected error (in perm(1)) to the lowest error (in
      %         perm(end)

      % sort provided y-values
      [ySort, yInd] = sort(yPredict);
      % sort provided variances, too
      sd2Sort = sd2Predict(yInd);
      % initializations
      n = length(ySort);
      expectedError = zeros(1,n);
      probMatrix = zeros(n,n);

      % calculate expected rank difference for each point
      for i = 1:n
        y = ySort(i);
        sd2 = sd2Sort(i);
        probs = zeros(1,n);

        % probabilities of reaching the values of the other points
        % prob(2:end) = normcdf(ySort([1:(i-1) (i+1):n]), y, sd2);
        probY = 1-normcdf(ySort, y, sd2);
        % probability of being first (unless previously first --
        % something is then added to this value)
        probs(1) = 1-probY(1);
        % probY from the last iteration of the next cycle
        lastProb = probY(1);

        % probabilities of other positions
        position = 1;
        for j = 1:n
          probs(position) = probs(position) + lastProb - probY(j);
          lastProb = probY(j);
          if (i ~= j)
            % find out what is the current inspected permutation
            permutation = 1:n;
            permutation(i) = [];
            permutation = [permutation(1:j-1) i permutation(j:end)];
            % save this contribution to the final expected error
            % of this (i)-th point
            expectedError(i) = expectedError(i) + probs(position) * obj.rankFunc(permutation, mu);

            % we will move on to calculation of the next probability
            position = position + 1;
          end
        end
        % fill also the last probability
        probs(position) = probs(position) + probY(end);
        % save the probabilities into the final probability matrix
        probMatrix(i,:) = probs;

        % Debug:
        % fprintf('Expected error of [%d] is %f.\n', i, expectedError(i));
      end

      % Debug:
      fprintf('Ranking of errors of sorted points: %s\n', num2str(ranking(-expectedError)));
      % fprintf('Final permutation of original points: %s\n', num2str(perm'));

      % return the final permutation of the points from the highest expected error
      % to the lowest (according to (-1)*expectedError sorted according to
      % inverse sort defined by yInd, which is eqal to ranking(yPredict)
      % yRnk = ranking(yPredict);
      % eRnk = ranking(-expectedError);

      [~, eInd] = sort(-expectedError);
      % the final order of points to reevaluate is following
      perm = yInd(eInd);

      % Debug:
      yRnk = ranking(yPredict);
      assert(all(yRnk(perm)' == eInd), 'Ranking of ''perm'' and ''expectedError'' is not same!');
    end


    function nInit = getProbNOrigInit(obj)
      fracN = obj.nOrigInit - floor(obj.nOrigInit);
      plus = 0;
      if (fracN > 0)
        plus = (rand() < fracN);
      end
      nInit = floor(obj.nOrigInit) + plus;
    end

    function succ = trainModelOrUseLast(obj, xTrain, yTrain, cmaesState, sampleOpts)
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
