classdef DoubleTrainedEC < EvolutionControl & Observable
  properties 
    model
    pop
    cmaesState
    counteval
    
    origRatioUpdater
    restrictedParam
    useDoubleTraining
    retrainedModel
    stats
    archive
    nPresampledPoints
    surrogateOpts
  end
  
  methods 
    function obj = DoubleTrainedEC(surrogateOpts, varargin)
    % constructor
      obj@Observable();
      obj.model = [];

      surrogateOpts.updaterType = defopts(surrogateOpts, 'updaterType', 'none');
      surrogateOpts.updaterParams = defopts(surrogateOpts, 'updaterParams', {});
      obj.origRatioUpdater = OrigRatioUpdaterFactory.createUpdater(obj, surrogateOpts);
      obj.restrictedParam = defopts(surrogateOpts, 'evoControlRestrictedParam', 0.1);

      obj.useDoubleTraining = defopts(surrogateOpts, 'evoControlUseDoubleTraining', true);
      obj.pop = [];
      obj.surrogateOpts = surrogateOpts;
      obj.stats = struct( ...
          'fmin', NaN, ...              % minimal original fitness in population
          'rmseReeval', NaN, ...        % RMSE of the re-evaluated point(s)
          'kendallReeval', NaN, ...     % Kendall's corr. of the re-evaluated point(s)
          'rankErrReeval', NaN, ...     % rank error of popul. with re-evaluated point(s)
          'rankErr2Models', NaN, ...    % rank error between prediction of two models
          'rmseValid', NaN, ...         % RMSE of the (2nd) model on the validation set
          'kendallValid', NaN, ...      % Kendall of the (2nd) model on the validation set
          'rankErrValid', NaN, ...      % rank error between true fitness and model pred.
          ...                           %       on the validation set
          'lastUsedOrigRatio', NaN, ... % restricted param which was used (last) in the last generation
          'adaptRankDiff', NaN, ...     % last measured rankDiff during update()
          'adaptGain', NaN, ...         % gain of original ratio (to be converted via min/max)
          'adaptNewRatio', NaN ...      % new value of ratio (to be history-weighted)
          );
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control
      
      % initialization
      lambda = cmaesState.lambda;       % this is needed due to myeval
      dim    = cmaesState.dim;          % this is needed due to myeval
      obj.cmaesState = cmaesState;
      obj.archive    = archive;
      obj.counteval  = counteval;
      obj.retrainedModel = [];
      obj.stats.nDataInRange = NaN;
      
      % prepare the final population to be returned to CMA-ES
      obj.pop = Population(lambda, dim);
      
      newModel = ModelFactory.createModel(obj.surrogateOpts.modelType, obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');

      if (isempty(newModel))
        % model could not be created :(. Use the standard CMA-ES.
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
            = obj.finalizeGeneration(sampleOpts, varargin);
        return;
      end
      
      minTrainSize = newModel.getNTrainData();

      nArchivePoints = myeval(obj.surrogateOpts.evoControlTrainNArchivePoints);
      [xTrain, yTrain, nData] = obj.archive.getDataNearPoint(nArchivePoints, ...
          obj.cmaesState.xmean', obj.surrogateOpts.evoControlTrainRange, ...
          obj.cmaesState.sigma, obj.cmaesState.BD);
      obj.stats.nDataInRange = nData;
      
      % Do pre-sample
      [ok, y, arx, x, arz, ~, obj.counteval, xTrain, yTrain] = ...
          presample(minTrainSize, obj.cmaesState, obj.surrogateOpts, sampleOpts, ...
          obj.archive, obj.counteval, xTrain, yTrain, varargin{:});
      obj.nPresampledPoints = size(x, 2);
      phase = 0;        % pre-sampled points
      obj.pop = obj.pop.addPoints(x, y, arx, arz, obj.nPresampledPoints, phase);

      if (~ok)
        % not enough data for training model ==> use original pop
        % TODO: try the old model instead just orig-evaluating all the population
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
            = obj.finalizeGeneration(sampleOpts, varargin);
        return;
      end

      % train the model 
      newModel = newModel.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
      if (~newModel.isTrained())
        % model cannot be trained :( -- return with orig-evaluated population
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
            = obj.finalizeGeneration(sampleOpts, varargin);
        % TODO: try the old model instead just orig-evaluating all the population
        return;
      else
        obj.model = newModel;
      end

      nLambdaRest = lambda - obj.nPresampledPoints;

      % sample new points
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, nLambdaRest, obj.cmaesState, sampleOpts);
      [modelOutput, yExtendModel] = obj.model.getModelOutput(xExtendValid');

      isEvaled = false(1, nLambdaRest);
      yOrig    = NaN(1, nLambdaRest);
      notEverythingEvaluated = true;

      while (notEverythingEvaluated)

        nPoints = ceil(nLambdaRest * obj.restrictedParam) - sum(isEvaled);
        obj.stats.lastUsedOrigRatio = obj.restrictedParam;
        % Debug:
        % fprintf('ratio: %.2f | nPoints: %d | iter: %d\n', obj.restrictedParam, nPoints, obj.cmaesState.countiter);

        reevalID = false(1, nLambdaRest);
        reevalID(~isEvaled) = obj.choosePointsForReevaluation(nPoints, ...
            xExtend(:, ~isEvaled), modelOutput(~isEvaled), yExtendModel(~isEvaled));
        xToReeval = xExtendValid(:, reevalID);
        nToReeval = sum(reevalID);

        % original-evaluate the chosen points
        [yNew, xNew, xNewValid, zNew, obj.counteval] = ...
            sampleCmaesOnlyFitness(xExtend(:, reevalID), xToReeval, zExtend(:, reevalID), ...
            obj.cmaesState.sigma, nToReeval, obj.counteval, obj.cmaesState, sampleOpts, ...
            varargin{:});
        xExtendValid(:, reevalID) = xNewValid;
        xExtend(:, reevalID) = xNew;
        zExtend(:, reevalID) = zNew;
        yOrig(reevalID) = yNew;
        isEvaled = isEvaled | reevalID;
        % Debug:
        % fprintf('counteval: %d\n', obj.counteval)

        phase = 1;        % re-evaluated points
        obj.pop = obj.pop.addPoints(xNewValid, yNew, xNew, zNew, nToReeval, phase);

        % update the Archive
        obj.archive.save(xNewValid', yNew', obj.cmaesState.countiter);

        % re-train the model again with the new original-evaluated points
        if ~all(isEvaled)
          xTrain = [xTrain; xNewValid'];
          yTrain = [yTrain; yNew'];
          obj.retrainedModel = obj.model.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
          if (obj.useDoubleTraining && obj.retrainedModel.isTrained())
            % origRatio adaptivity
            if (~isempty(obj.origRatioUpdater.lastUpdateGeneration) ...
                && obj.origRatioUpdater.lastUpdateGeneration > obj.cmaesState.countiter)
              % internal CMA-ES restart, create a new origRatioUpdater
              obj.origRatioUpdater = OrigRatioUpdaterFactory.createUpdater(obj, obj.surrogateOpts);
            end
            yFirstModel  = obj.model.predict(xExtendValid');
            yExtendModel = obj.retrainedModel.predict(xExtendValid');
            obj.restrictedParam = obj.origRatioUpdater.update(...
                yFirstModel', yExtendModel', dim, lambda, obj.cmaesState.countiter, obj);
            obj.stats.adaptRankDiff = obj.origRatioUpdater.rankDiffs(end);
            obj.stats.adaptGain = obj.origRatioUpdater.gain;
            obj.stats.adaptNewRatio = obj.origRatioUpdater.newRatio;
            % Debug:
            % fprintf('OrigRatio: %f\n', obj.origRatioUpdater.getLastRatio(obj.cmaesState.countiter));

            % Debug:
            % fprintf('UPDATE: ratio: %.2f | rankDiff: %.2f | iter: %d\n', obj.restrictedParam, obj.origRatioUpdater.rankDiffs(end), obj.cmaesState.countiter);
          else
            % Debug:
            % fprintf('DoubleTrainedEC: The new model could (is not set to) be trained, using the not-retrained model.\n');
          end
        end
      
        notEverythingEvaluated = (floor(lambda * obj.restrictedParam) > sum(isEvaled));
      end

      if (~all(isEvaled))
        phase = 2;      % model-evaluated rest of the population
        obj.pop = obj.pop.addPoints(xExtendValid(:, ~isEvaled), yExtendModel(~isEvaled), ...
            xExtend(:, ~isEvaled), zExtend(:, ~isEvaled), 0, phase);
      end
      
      assert(obj.pop.nPoints == lambda, 'There are not yet all lambda points prepared, but they should be!');

      % sort the returned solutions (the best to be first)
      [obj.pop, ~] = obj.pop.sort;

      % save the resulting re-evaluated population as the returning parameters
      [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
          = obj.finalizeGeneration(sampleOpts, varargin);
    end


    function [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = finalizeGeneration(obj, sampleOpts, varargin)
      % fill the rest of the population with original evaluations
      [obj, fitness_raw, arx, arxvalid, arz, counteval] = ...
          obj.fillPopWithOrigFitness(sampleOpts, varargin);
      origEvaled = obj.pop.origEvaled;

      % calculate statistics
      obj.stats.rmseReeval     = NaN; % RMSE of the re-evaluated point(s)
      obj.stats.kendallReeval  = NaN; % Kendall's corr. of the re-evaluated point(s)
      obj.stats.rankErrReeval  = NaN; % rank error of popul. with re-evaluated point(s)
      obj.stats.rankErr2Models = NaN; % rank error between prediction of two models
      obj.stats.rmseValid      = NaN; % RMSE of the (2nd) model on the validation set
      obj.stats.kendallValid   = NaN; % Kendall of the (2nd) model on the validation set
      obj.stats.rankErrValid   = NaN; % rank error between true fitness and model pred.
      obj.stats.fmin = min(obj.pop.y(1,(obj.pop.origEvaled == true)));

      % model-related statistics
      if (~isempty(obj.model) && obj.model.isTrained() ....
          && ~all(obj.pop.origEvaled))
        % predict the population by the first model
        yModel1 = obj.model.predict(obj.pop.x');

        % calculate RMSE, Kendall's coeff. and ranking error
        % between the original fitness and the first model's values
        % of the re-evaluated point(s), i.e. (phase == 1)
        [obj.stats.rmseReeval, obj.stats.kendallReeval, obj.stats.rankErrReeval] ...
            = obj.reevalStatistics(yModel1);

        % get ranking error between the first and the second model
        % (if the second model is trained)
        obj.stats.rankErr2Models = obj.retrainStatistics(yModel1);

        % independent validation set statistics
        [obj.stats.rmseValid, obj.stats.kendallValid, obj.stats.rankErrValid, lastModel] ...
            = obj.validationStatistics(sampleOpts);

        % shift the f-values:
        %   if the model predictions are better than the best original value
        %   in the model's dataset, shift ALL (!) function values
        %   Note: - all values have to be shifted in order to preserve predicted
        %           ordering of values
        %         - small constant is added because of the rounding errors
        %           when numbers of different orders of magnitude are summed
        fminDataset = min(lastModel.dataset.y);
        fminModel = obj.pop.getMinModeled;
        diff = max(fminDataset - fminModel, 0);
        obj.pop = obj.pop.shiftY(1.000001*diff);
        fitness_raw = obj.pop.y;
      end

      obj.notify_observers();

      % for backwards compatibility
      [minstd minstdidx] = min(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      [maxstd maxstdidx] = max(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      surrogateStats = [obj.stats.rmseValid, obj.stats.kendallValid, obj.cmaesState.sigma, ...
        max(obj.cmaesState.diagD)/min(obj.cmaesState.diagD), minstd, maxstd, ...
        obj.stats.rankErr2Models];
    end
    

    function [obj, yNew, xNew, xNewValid, zNew, counteval] = fillPopWithOrigFitness(obj, sampleOpts, varargin)
      %Fill the rest of the current population 'pop' (of class Population) with 
      % the original-evaluated individuals
      nToEval = obj.cmaesState.lambda - sum(obj.pop.nPoints);
      if (nToEval > 0)
        [y, arx, x, arz, obj.counteval] = sampleCmaes(obj.cmaesState, sampleOpts, ...
            nToEval, obj.counteval, varargin{:});
        obj.archive.save(x', y', obj.cmaesState.countiter);

        phase = 3;      % original-evaluated rest of the population
        obj.pop = obj.pop.addPoints(x, y, arx, arz, nToEval, phase);
      end
      obj.pop.sort();
      yNew = obj.pop.y;
      xNewValid = obj.pop.x;
      xNew = obj.pop.arx;
      zNew = obj.pop.arz;
      counteval = obj.counteval;
    end

    function reevalID = choosePointsForReevaluation(obj, nPoints, xExtend, modelOutput, yExtendModel)
    % choose points with low confidence to reevaluate
    %
    % returns:
    %   reevalID  --  bool vector which points from xExtend to reevaluate
      if any(strcmpi(obj.model.predictionType, {'sd2', 'poi', 'ei'}))
        % higher criterion is better (sd2, poi, ei)
        [~, pointID] = sort(modelOutput, 'descend');
      elseif (strcmpi(obj.model.predictionType, 'expectedrank'))
        % TODO it should work (yExtendModel, modelOutput) instead of re-predict the points again
        [y_m, sd2_m] = obj.model.predict(xExtend');
        pointID = expectedRankDiff(y_m, sd2_m, obj.cmaesState.mu, @errRankMuOnly);
        % Debug:
        % y_r = ranking(y_m);
        % fprintf('  Expected permutation of sorted f-values: %s\n', num2str(y_r(pointID)'));
      else
        % lower criterion is better (fvalues, lcb, fpoi, fei)
        [~, pointID] = sort(yExtendModel, 'ascend');
      end
      nLambdaRest = size(xExtend, 2);
      reevalID = false(1, nLambdaRest);
      assert(obj.origRatioUpdater.getLastRatio(obj.cmaesState.countiter) >= 0 && obj.origRatioUpdater.getLastRatio(obj.cmaesState.countiter) <= 1, 'origRatio out of bounds [0,1]');
      reevalID(pointID(1:nPoints)) = true;
    end


    function [rmse, kendall, rankErr] = reevalStatistics(obj, yModel1)
      % calculate RMSE and possibly Kendall's coeff. of the re-evaluated point(s)
      % (phase == 1)
      phase1 = (obj.pop.phase == 1);
      yReeval = obj.pop.y(1,phase1);
      yReevalModel = yModel1(phase1);
      rmse = sqrt(sum((yReevalModel' - yReeval).^2))/length(yReeval);
      kendall = corr(yReevalModel, yReeval', 'type', 'Kendall');

      [~, sort1] = sort(yModel1);
      yModel2 = yModel1;
      yModel2(phase1) = yReeval;
      ranking2   = ranking(yModel2);
      rankErr = errRankMuOnly(ranking2(sort1), obj.cmaesState.mu);

      % Debug:
      % fprintf('  model: %d pSmpls, reeval %d pts, RMSE= %.2e, Kendl= %.2f, rankErr= %.3f\n', ...
      %     obj.nPresampledPoints, length(yReeval), rmse, kendall, rankErr);
    end


    function rankErr = retrainStatistics(obj, yModel1)
    % get ranking error between the first and the second model
    % (if the second model is trained)
      if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
        yModel1AfterPresample = yModel1(obj.pop.phase ~= 0);
        xAfterPresample = obj.pop.x(:,obj.pop.phase ~= 0);
        yModel2AfterPresample = obj.retrainedModel.predict(xAfterPresample');
        [~, sort1] = sort(yModel1AfterPresample);
        ranking2   = ranking(yModel2AfterPresample);
        rankErr = errRankMuOnly(ranking2(sort1), obj.cmaesState.mu);
        % Debug:
        % fprintf('  2 models rank error: %.3f                  %s\n', rankErr, decorateKendall(1-rankErr*2));
      else
        % ranking error between two models cannot be calculated :(
        rankErr = NaN;
      end
    end


    function [rmse, kendall, errRank, lastModel] = validationStatistics(obj, sampleOpts)
    % generate independent validation population of points,
    % evaluate them with the original BBOB function (if we know it)
    % and calculate statistics on this independent set
    %
    % lastModel -- return the last valid model (either first or second)

      rmse = NaN; kendall = NaN; errRank = NaN;
      % do we have access to the original BBOB fitness?
      if (~ isfield(obj.surrogateOpts.modelOpts, 'bbob_func'))
        return;
      end
      if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
        % statistics according to the retrained model
        lastModel = obj.retrainedModel;
      elseif (~isempty(obj.model) && obj.model.isTrained())
        % statistics according to the first model
        lastModel = obj.model;
      else
        % we do not have any valid model
        return;
      end

      [~, xValidTest, ~] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, obj.cmaesState.lambda, obj.cmaesState, sampleOpts);
      preciseModel = ModelFactory.createModel('bbob', obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');
      yTest = preciseModel.predict(xValidTest');
      yPredict = lastModel.predict(xValidTest');

      kendall = corr(yPredict, yTest, 'type', 'Kendall');
      rmse = sqrt(sum((yPredict - yTest).^2))/length(yPredict);

      [~, sort1] = sort(yTest);
      ranking2   = ranking(yPredict);
      errRank = errRankMuOnly(ranking2(sort1), obj.cmaesState.mu);

      % Debug:
      % fprintf('  test RMSE= %.2e, Kendall= %.3f, rankErr= %.3f %s\n', ...
      %     rmse, kendall, errRank, decorateKendall(kendall));
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
