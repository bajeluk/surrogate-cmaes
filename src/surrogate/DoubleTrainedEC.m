classdef DoubleTrainedEC < EvolutionControl & Observable
  properties 
    model
    pop
    cmaesState
    counteval
    
    restrictedParam
    useDoubleTraining
    retrainedModel
    stats
    archive
    nPresampledPoints
    surrogateOpts
    newModel
    modelArchive
    modelArchiveGenerations
    modelArchiveLength
    acceptedModelAge            % how many generations old model is still OK (0 == only current)
    modelAge                    % age of model in the number of generations (0 == current model)
    oldModelAgeForStatistics    % age of model for gathering statistics of old models
  end
  
  methods 
    function obj = DoubleTrainedEC(surrogateOpts, varargin)
    % constructor
      obj@Observable();
      obj.model = [];
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
          'rankErrValid', NaN, ...      % rank error between true fitn. and model pred. on the validation set
          'rmseOldModel', NaN, ...      % RMSE of old model on future orig points
          'kendallOldModel', NaN, ...   % Kendall of old model on future orig points
          'rankErrOldModel', NaN, ...   % rank error between true fitness and model pred.
          'ageOldModel', -1 ...         % age of the old model which used for statistics
          );
      obj.modelAge = 0;

      % fixed settings:
      obj.oldModelAgeForStatistics = [3:5];
      obj.modelArchiveLength = 5;

      % other initializations:
      obj.modelArchive = cell(1, obj.modelArchiveLength);
      obj.modelArchiveGenerations = (-1)*ones(1, obj.modelArchiveLength);
      obj.acceptedModelAge = defopts(surrogateOpts, 'evoControlAcceptedModelAge', 2);
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
      obj.modelAge = 0;
      
      % prepare the final population to be returned to CMA-ES
      obj.pop = Population(lambda, dim);
      
      obj.newModel = ModelFactory.createModel(obj.surrogateOpts.modelType, obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');

      if (isempty(obj.newModel))
        [obj, ok] = tryOldModel();
        if (~ok)
          % model could not be created nor older is usable :(. Use the standard CMA-ES.
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
              = obj.finalizeGeneration(sampleOpts, varargin);
          return;
        end
      end

      % if a new model is used, find appropriate training set and train it:
      if (obj.modelAge == 0)

      minTrainSize = obj.newModel.getNTrainData();

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
        % not enough data for training model
        [obj, ok] = obj.tryOldModel();
        if (~ok)
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
              = obj.finalizeGeneration(sampleOpts, varargin);
          return;
        end
      end

      % train the model 
      obj.newModel = obj.newModel.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
      if (obj.newModel.isTrained())
        obj = obj.updateModelArchive(obj.newModel, obj.modelAge);
      else
        [obj, ok] = tryOldModel();
        if (~ok)
          % model cannot be trained :( -- return with orig-evaluated population
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
              = obj.finalizeGeneration(sampleOpts, varargin);
          return;
        end
      end

      end  % if (obj.modelAge == 0)

      obj.model = obj.newModel;

      nLambdaRest = lambda - obj.nPresampledPoints;

      % sample new points
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, nLambdaRest, obj.cmaesState, sampleOpts);
      [modelOutput, yExtendModel] = obj.model.getModelOutput(xExtend');

      reevalID = obj.choosePointsForReevaluation(xExtend, modelOutput, yExtendModel);
      xToReeval = xExtendValid(:, reevalID);
      nToReeval = sum(reevalID);

      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, obj.counteval] = ...
          sampleCmaesOnlyFitness(xExtend(:, reevalID), xToReeval, zExtend(:, reevalID), ...
          obj.cmaesState.sigma, nToReeval, obj.counteval, obj.cmaesState, sampleOpts, ...
          varargin{:});
      % Debug:
      % fprintf('counteval: %d\n', obj.counteval)

      phase = 1;        % re-evaluated points
      obj.pop = obj.pop.addPoints(xNewValid, yNew, xNew, zNew, nToReeval, phase);

      % update the Archive
      obj.archive.save(xNewValid', yNew', obj.cmaesState.countiter);

      % re-train the model again with the new original-evaluated points
      if ~all(reevalID)
        xTrain = [xTrain; xNewValid'];
        yTrain = [yTrain; yNew'];
        obj.retrainedModel = obj.model.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
        if (obj.useDoubleTraining && obj.retrainedModel.isTrained())
          obj = obj.updateModelArchive(obj.retrainedModel, obj.modelAge);
          yReeval = obj.retrainedModel.predict((xExtendValid(:, ~reevalID))');
        else
          % Debug:
          % fprintf('DoubleTrainedEC: The new model could (is not set to) be trained, using the not-retrained model.\n');
          yReeval = yExtendModel(~reevalID);
        end

        phase = 2;      % model-evaluated rest of the population
        obj.pop = obj.pop.addPoints(xExtendValid(:, ~reevalID), yReeval, ...
            xExtend(:, ~reevalID), zExtend(:, ~reevalID), 0, phase);
      end
              
      assert(obj.pop.nPoints == lambda, 'There are not yet all lambda points prepared, but they should be!');

      % sort the returned solutions (the best to be first)
      [obj.pop, ~] = obj.pop.sort;

      % save the resulting re-evaluated population as the returning parameters
      [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
          = obj.finalizeGeneration(sampleOpts, varargin);
    end


    function obj = updateModelArchive(obj, newModel, modelAge)
      % update the modelArchive with the current new model
      countiter = obj.cmaesState.countiter;

      if (obj.modelArchiveGenerations(1) < (countiter - modelAge))
        % there's an old model in the first position ==> shift old
        % models to the history
        obj.modelArchive(2:end) = obj.modelArchive(1:(end-1));
        obj.modelArchiveGenerations(2:end) = obj.modelArchiveGenerations(1:(end-1));
        % clear the first position
        obj.modelArchive{1} = [];
        obj.modelArchiveGenerations(1) = -1;
      end

      if (newModel.isTrained())
        % the obj.newModel should be usable and newer than what we have
        % in the first position, so save it there
        obj.modelArchive{1} = newModel;
        obj.modelArchiveGenerations(1) = (countiter - modelAge);
      end
    end


    function [obj, ok] = tryOldModel(obj)
      % try to load an appropriate old model
      % if successful, load it as the newModel
      [oldModel, age] = obj.getOldModel(0:obj.acceptedModelAge);
      if (~isempty(oldModel))
        obj.newModel = oldModel;
        obj.modelAge = age;
        ok = true;
      else
        ok = false;
      end
    end


    function [oldModel, modelAge] = getOldModel(obj, generationDiffs)
      % return an old model from modelArchive from the generations
      %   countiter - [generationDiffs]
      % Returns the first such model found, or [] if none is found
      %
      % Examples:
      % m = dec.getOldModel(3:5)        % returns the youngest model 3--5 generations old
      % m = dec.getOldModel(0)          % returns the model from current generation

      countiter = obj.cmaesState.countiter;
      modelGenerations = countiter - generationDiffs;
      modelGenerations = modelGenerations(modelGenerations > 0);

      oldModel = []; modelAge = -1;
      if (isempty(modelGenerations))
        % no sensible generations given, return []
        return;
      end

      % go through modelArchive and look whether the models (in
      % historical order) are not from the specified generations
      for i = 1:obj.modelArchiveLength
        ind = find(modelGenerations == obj.modelArchiveGenerations(i), 1, 'first');
        if (ind)
          oldModel = obj.modelArchive{i};
          modelAge = (countiter - obj.modelArchiveGenerations(i));
          return;
        end
      end
      % no model from specified generations found, return []
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
      obj.stats.rmseOldModel   = NaN; % RMSE of old model on future orig points
      obj.stats.kendallOldModel = NaN; % Kendall of old model on future orig points
      obj.stats.rankErrOldModel = NaN; % rank error between true fitness and model pred.
      obj.stats.ageOldModel    = -1;
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

        [obj.stats.rmseOldModel, obj.stats.kendallOldModel, obj.stats.rankErrOldModel, obj.stats.ageOldModel] ...
            = obj.oldModelStatistics();

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

    function reevalID = choosePointsForReevaluation(obj, xExtend, modelOutput, yExtendModel)
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
        [~, pointID] = sort(modelOutput, 'ascend');
      end
      nLambdaRest = size(xExtend, 2);
      reevalID = false(1, nLambdaRest);
      assert(obj.restrictedParam >= 0 && obj.restrictedParam <= 1, 'evoControlRestrictedParam out of bounds [0,1]');
      nToReeval = ceil(nLambdaRest * obj.restrictedParam);
      reevalID(pointID(1:nToReeval)) = true;
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


    function [rmse, kendall, rankErr, age] = oldModelStatistics(obj)
      % return statistics of a several-generations-old model measured
      % on new data from following generations (generations after the
      % model was trained)
      rmse = NaN; kendall = NaN; rankErr = NaN;

      [oldModel, age] = obj.getOldModel(obj.oldModelAgeForStatistics);
      if (~isempty(oldModel))
        countiter = obj.cmaesState.countiter;
        oldModelGeneration = oldModel.trainGeneration;
        testGenerations = [(oldModel.trainGeneration+1):countiter];
        [xOrig, yOrig] = obj.archive.getDataFromGenerations(testGenerations);
        if (~isempty(yOrig))
          yPredict = oldModel.predict(xOrig);
          rmse = sqrt(sum((yPredict - yOrig).^2))/length(yPredict);
          kendall = corr(yPredict, yOrig, 'type', 'Kendall');
          [~, sort1] = sort(yPredict);
          ranking2   = ranking(yOrig);
          rankErr = errRankMuOnly(ranking2(sort1), obj.cmaesState.mu);
        end
      end
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
