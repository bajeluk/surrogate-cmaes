classdef LinQuadEC < EvolutionControl & Observable
%
% TODO:
% [ ] remove updaterParams and use DTAdaptive_* parameters instead
% [ ] rename 'restrictedParam' to 'origRatio'
% [ ] when preselection, remove not the first points, but be a bit more clever (e.g. the worst predicted points...?)
%
  properties
    model
    pop
    cmaesState
    counteval

    origRatioUpdater
    restrictedParam
    useDoubleTraining
    maxDoubleTrainIterations
    minPointsForExpectedRank    % minimal # of points in generation for expectedRank criterion to be allowed to use, otherwise sd2 criterion is used
    retrainedModel
    stats
    usedUpdaterState            % Updater's state variables updated in the last generation
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
    isTrainSuccess
    trainRange
    origPointsRoundFcn          % function computing number of original-evaluated points from origRatio
    nBestPoints                 % the number of points with the best predicted f-value to take every generation
    usedBestPoints              % how many best-predicted points was really orig-evaluated
    preselectionPopRatio        % how many times larger population should be used for preselection
    validationGenerationPeriod  % the number of generations between "validation generations" + 1, validation generation is used when 'mod(g, validationGenerationPeriod) == 0'; see validationPopSize
    validationPopSize           % the minimal number of points to be orig-evaluated in validation generation
    
    individualToInject
    useInject
  end

  methods
    function obj = LinQuadEC(surrogateOpts, varargin)
    % constructor
      obj@Observable();
      obj.pop = [];
      obj.counteval = 0;
      obj.surrogateOpts = surrogateOpts;

      % Model Archive fixed settings and initialization
      obj.modelArchiveLength = defopts(surrogateOpts, 'evoControlModelArchiveLength', 5);
      obj.oldModelAgeForStatistics = defopts(surrogateOpts, 'evoControlOldModelAgeForStatistics', ...
          3:min(5, obj.modelArchiveLength));
      obj.modelArchive = cell(1, obj.modelArchiveLength);
      obj.modelArchiveGenerations = nan(1, obj.modelArchiveLength);
      obj.modelAge = 0;
      obj.isTrainSuccess = false;
      
      obj.acceptedModelAge = defopts(surrogateOpts, 'evoControlAcceptedModelAge', 2);
      obj.origPointsRoundFcn = str2func(defopts(surrogateOpts, 'evoControlOrigPointsRoundFcn', 'ceil'));
      obj.trainRange = defopts(surrogateOpts, 'evoControlTrainRange', 10);
      
      obj.useInject = defopts(surrogateOpts, 'evoControlUseInject', true);

      % statistics
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
          'normKendallOldModel', NaN, ... % Kendall of old model normed to [0,1]
          'ageOldModel', NaN, ...       % age of the old model which used for statistics
          'nDataOldModel', 0, ...       % the number of data points from archive for old model statistics
          'lastUsedOrigRatio', NaN, ... % restricted param which was used (last) in the last generation
          'adaptErr', NaN, ...     % last measured rankDiff during update()
          'adaptGain', NaN, ...         % gain of original ratio (to be converted via min/max)
          'adaptSmoothedErr', NaN ...   % smoothed error value used before fed into transfer function
          );
      obj.usedUpdaterState = struct( ...
          'gain', NaN, ...
          'err', NaN, ...
          'smoothedErr', NaN ...
          );
      
      obj.individualToInject = NaN
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
      obj.isTrainSuccess = false;
      obj.usedBestPoints = 0;
      if (ischar(obj.trainRange))
        obj.trainRange = round(myeval(obj.trainRange));
      end

      % prepare the final population to be returned to CMA-ES
      obj.pop = Population(lambda, dim);      

      obj.newModel = ModelFactory.createModel(obj.surrogateOpts.modelType, obj.surrogateOpts.modelOpts, obj.cmaesState.xmean', obj.model);

      if (isempty(obj.newModel))
        [obj, ok] = obj.tryOldModel();
        if (~ok)
          % model could not be created nor older is usable :(. Use the standard CMA-ES.
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
              = obj.finalizeGeneration(sampleOpts, varargin{:});
          return;
        end
      end

      % if a new model is used, find appropriate training set and train it:
      if (isempty(obj.model))

        minTrainSize = obj.newModel.getNTrainData();

        if isfield(obj.surrogateOpts.modelOpts, 'trainsetSizeMax')
            nArchivePoints = myeval(obj.surrogateOpts.modelOpts.trainsetSizeMax);
        else
            nArchivePoints = myeval(obj.surrogateOpts.evoControlTrainNArchivePoints);
        end
        
        
        % Choose data from Archive -- the points from this 'getDataNearPoint() 
        % are used only iff modelOpts.trainsetType == 'parameters', otherwise the points
        % are selected later in Model.train() -> ... -> Archive.getTrainsetData()
        [xTrain, yTrain, nData] = obj.archive.getDataNearPoint(nArchivePoints, ...
            obj.cmaesState.xmean', obj.trainRange, ...
            obj.cmaesState.sigma, obj.cmaesState.BD);
        obj.stats.nDataInRange = nData;

        % Do pre-sample
        [ok, y, arx, x, arz, ~, obj.counteval] = ...
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
                = obj.finalizeGeneration(sampleOpts, varargin{:});
            return;
          end
        end
      end

     
      % sample new population of lambda points out of which will be chosen
      % later in this generation
      nLambdaRest = lambda - obj.pop.nPoints;
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, nLambdaRest, obj.cmaesState, sampleOpts, obj.individualToInject);
      
      evaluated = obj.pop.nPoints;
      k = evaluated;
      
      
      % add these points into Population without valid f-value (marked as notEvaluated)
      phase = 4;
      obj.pop = obj.pop.addPoints(xExtendValid, NaN(1,nLambdaRest), xExtend, zExtend, 0, phase);

      if (isempty(obj.model))
        % (first) model training
        X_tr = [xTrain; obj.pop.getOriginalX()'];
        y_tr = [yTrain; obj.pop.getOriginalY()'];
        obj.newModel = obj.newModel.train(X_tr, y_tr, obj.cmaesState, sampleOpts, obj.archive, obj.pop);
        if (obj.newModel.isTrained())
          obj = obj.updateModelArchive(obj.newModel, obj.modelAge);
          obj.model = obj.newModel;
        else
          [obj, ok] = obj.tryOldModel();
          if (~ok)
            % model cannot be trained :( -- return with orig-evaluated population
            [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
                = obj.finalizeGeneration(sampleOpts, varargin{:});
            return;
          end
        end
      end  % if (obj.modelAge == 0)
      
      
      %modelPredictions = obj.newModel.modelPredict(obj.pop.x');
      %[sortedModelPredictions, sortedModelPredictionsIndexes] = sort(modelPredictions);
      
      if (k == 0) 
        k = floor(1 + lambda*0.02);
      end

      while(evaluated < lambda)
          [sortedModelPredictions, sortedModelPredictionsIndexes] = sort(obj.model.predict(obj.pop.x'));
          firstKIndexes = sortedModelPredictionsIndexes(1:k);
          i = 1;
          while evaluated < k
              idxToEval = firstKIndexes(i);
              i = i + 1;
              if ~obj.pop.isEvaled(idxToEval)
                  toEval = obj.pop.x(:, idxToEval);
                  [y, arx, x, arz, obj.counteval] = ...
                    sampleCmaesOnlyFitness(toEval, toEval, toEval, obj.cmaesState.sigma, 1, ...
                    obj.counteval, obj.cmaesState, sampleOpts, 'Archive', obj.archive, varargin{:});
                  obj.archive.save(x', y', obj.cmaesState.countiter);
                  evaluated = evaluated + 1;
                  
                  logicalIndexes = false(1, lambda);
                  logicalIndexes(idxToEval) = true;
                  obj.pop = obj.pop.updateYValue(x, y, 1, phase, logicalIndexes);
              end             
          end
          obj.archive.sortLast(k);
          
          obj.newModel = obj.newModel.train([], [], obj.cmaesState, sampleOpts, obj.archive, obj.pop);
          
          if obj.newModel.isTrained()
            obj.model = obj.newModel;
            err = obj.calcKendallError(k);
            if err >= 0.85
                break;
            end
          end
          
          %currIndex = currIndex + evalCount;
          k = ceil(k*1.5);
          if k > lambda
            k = lambda;
          end
      end
      
      if obj.useInject
        modelMinimum = obj.model.minimumX(obj.archive);
        rescaling_factor = normrnd(0, 1, [1, dim]);
        rescaling_factor = rescaling_factor .^ 2;
        rescaling_factor = sum(rescaling_factor);
        rescaling_factor = rescaling_factor .^ 0.5;
        rescaling_factor = rescaling_factor / obj.mahalanobisNorm(modelMinimum);
        obj.individualToInject = modelMinimum * (rescaling_factor / obj.cmaesState.sigma);
      end
      
      
      if (evaluated >= lambda)
        %All points were evaluated -> return orig evaluated solution
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
                = obj.finalizeGeneration(sampleOpts, varargin{:});
            return;
      else
          % Model is good enough to evaluate whole population
          modelPredictions = obj.model.predict(obj.pop.x');
          offset = min(obj.pop.getOriginalY) - min(modelPredictions);
          modelPredictions = modelPredictions + offset + 0.000001*abs(offset);
          obj.pop = obj.pop.updateYValue([], modelPredictions, 0, phase, true(1, lambda));
          %obj.pop.origEvaled = false(1, lambda);
          [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] ...
                 = obj.finalizeGeneration(sampleOpts, varargin{:});
          return;
          
      end

    end
    
    function res = mahalanobisNorm(obj, x)
        D = diag(obj.cmaesState.diagD);
        B = obj.cmaesState.BD * inv(D);
        res = sum(((B * x') ./ obj.cmaesState.diagD) .^2 ) .^ 0.5;
        res = res / obj.cmaesState.sigma;
    end
    
    function err = calcKendallError(obj, k)
        valuesToTestCnt = round(max(15, min(1.2 * k, 0.75 * obj.cmaesState.lambda)));
        [archiveSize, ~] = size(obj.archive.y);
        valuesToTestCnt = min(archiveSize, valuesToTestCnt);
        
        X = obj.archive.X(archiveSize:-1:(archiveSize - valuesToTestCnt + 1), :);
        modelPrediction = obj.model.predict(X);
        y = obj.archive.y(archiveSize:-1:(archiveSize - valuesToTestCnt + 1));
        try
            err = corr(modelPrediction, y, 'type', 'Kendall');
        catch
            display(err);
        end
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
        obj.modelArchiveGenerations(1) = NaN;
        obj.isTrainSuccess = true;
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

      oldModel = []; modelAge = NaN;
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
          obj.fillPopWithOrigFitness(sampleOpts, varargin{:});
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
      obj.stats.normKendallOldModel = NaN; % Kendall transformed to [0,1] error-range
      obj.stats.ageOldModel    = NaN;
      obj.stats.nDataOldModel  = 0;
      obj.stats.fmin = min(obj.pop.y(1,(obj.pop.origEvaled == true)));

      [obj.stats.rmseOldModel, obj.stats.kendallOldModel, ...
          obj.stats.normKendallOldModel, obj.stats.ageOldModel, ...
          obj.stats.nDataOldModel] = obj.oldModelStatistics();

      % model-related statistics
      if (~isempty(obj.model) && obj.model.isTrained())
        % predict the population by the first model
        
        if (any(obj.pop.origEvaled))
          % calculate RMSE, Kendall's coeff. and ranking error
          % between the original fitness and the first model's values
          % of the re-evaluated point(s), i.e. (phase == 1)
          
          
          %[obj.stats.rmseReeval, obj.stats.kendallReeval, obj.stats.rankErrReeval] ...
          %    = obj.reevalStatistics();

          % get ranking error between the first and the second model
          % (if the second model is trained)
          % obj.stats.rankErr2Models = obj.retrainStatistics(yModel1);
        end

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
        fminModel = obj.pop.getMinModeled;
        if (~isempty(fminModel))
          fminDataset = min(lastModel.getDataset_y());
          diff = max(fminDataset - fminModel, 0);
          obj.pop = obj.pop.shiftY(1.000001*diff);
          fitness_raw = obj.pop.y;
        end
      end

      obj.notify_observers();

      % for backwards compatibility
      [minstd, minstdidx] = min(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      [maxstd, maxstdidx] = max(obj.cmaesState.sigma*sqrt(obj.cmaesState.diagC));
      surrogateStats = [obj.stats.rmseValid, obj.stats.kendallValid, obj.cmaesState.sigma, ...
        max(obj.cmaesState.diagD)/min(obj.cmaesState.diagD), minstd, maxstd, ...
        obj.stats.rankErr2Models];
    end


    function [obj, yNew, xNew, xNewValid, zNew, counteval] = fillPopWithOrigFitness(obj, sampleOpts, varargin)
      %Fill the rest of the current population 'pop' (of class Population) with
      % the original-evaluated individuals
      [x, arz] = obj.pop.getNotEvaledX();
      nToEval = obj.pop.lambda - sum(obj.pop.isEvaled);
      if (nToEval > 0)
        nSample = nToEval - size(x, 2);
        if (nSample > 0)
          [arxFill, xFill, zFill] = ...
            sampleCmaesNoFitness(obj.cmaesState.sigma, nSample, obj.cmaesState, sampleOpts);
          x   = [x, xFill];
          arz = [arz, zFill];
          phase = 4;    % not-evaluated points
          obj.pop = obj.pop.addPoints(xFill, NaN(1,size(xFill,2)), arxFill, zFill, 0, phase);
        end
        [y, arx, x, arz, obj.counteval] = ...
            sampleCmaesOnlyFitness(x, x, arz, obj.cmaesState.sigma, nToEval, ...
            obj.counteval, obj.cmaesState, sampleOpts, 'Archive', obj.archive, varargin{:});
        obj.archive.save(x', y', obj.cmaesState.countiter);

        phase = 3;      % original-evaluated rest of the population
        obj.pop = obj.pop.updateYValue(x, y, nToEval, phase);
      end

      obj.pop = obj.pop.sort();
      yNew = obj.pop.y;
      xNewValid = obj.pop.x;
      xNew = obj.pop.arx;
      zNew = obj.pop.arz;
      counteval = obj.counteval;

      assert(~any(isnan(yNew)), 'Assertion failed: fillPopWithOrigFitness is about to return some NaN''s');
    end

    function reevalID = choosePointsForReevaluation(obj, pop, thisModel, nPoints)
    % choose 'nPoints' points with the highest value of the criterion for
    % original reevaluation
    %
    % returns:
    %   reevalID  --  bool vector of points from xExtend to reevaluate
    %
    % TODO:
    % [ ] consider feeding all the population into expectedRankDiff()
      notOrigEvaledX = pop.getNotOrigEvaledX();

      if any(strcmpi(thisModel.predictionType, {'sd2', 'poi', 'ei'}))
        % higher criterion is better (sd2, poi, ei)
        modelOutput = thisModel.getModelOutput(notOrigEvaledX');
        [~, pointID] = sort(modelOutput, 'descend');

      elseif (strcmpi(thisModel.predictionType, 'expectedrank'))
        ok = true;
        if (isempty(thisModel) || ~thisModel.isTrained())
          warning('No valid model for calculating expectedRankDiff(). Using "sd2" criterion.');
          ok = false;
        end
        if (size(notOrigEvaledX, 2) < obj.minPointsForExpectedRank)
          fprintf(2, 'expectedRankDiff(): #pop=%d < %d: using sd2 criterion\n', size(notOrigEvaledX, 2), obj.minPointsForExpectedRank);
          ok = false;
        end
        if (ok)
          mu = ceil(obj.cmaesState.mu * (size(notOrigEvaledX ,2) / obj.cmaesState.lambda));
          [pointID, errs] = expectedRankDiff(thisModel, notOrigEvaledX, mu);
          if (~ sum(errs >= eps) > (size(notOrigEvaledX,2)/2))
            warning('exptectedRankDiff() returned more than lambda/2 points with zero expected rankDiff error. Using "sd2" criterion.');
            ok = false;
          end
        end
        if (~ok)
          [~, sd2] = thisModel.predict(notOrigEvaledX');
          [~, pointID] = sort(sd2, 'descend');
        end
        % Debug:
        % y_r = ranking(y_m);
        % fprintf('  Expected permutation of sorted f-values: %s\n', num2str(y_r(pointID)'));

      else
        % lower criterion is better (fvalues, lcb, fpoi, fei)
        modelOutput = thisModel.getModelOutput(notOrigEvaledX');
        [~, pointID] = sort(modelOutput, 'ascend');
      end

      notOrigEvaledID = find(~pop.origEvaled);
      reevalID = false(1, pop.lambda);
      reevalID(notOrigEvaledID(pointID(1:nPoints))) = true;
      %
      % Check the value of origRatio
      assert(obj.origRatioUpdater.getLastRatio() >= 0 && obj.origRatioUpdater.getLastRatio() <= 1, 'origRatio out of bounds [0,1]');
    end


    function [rmse, kendall, rankErr] = reevalStatistics(obj)
      % calculate RMSE and possibly Kendall's coeff. of the re-evaluated point(s)
      % (phase == 1)
      if (~any(obj.pop.origEvaled))
        rmse = NaN; kendall = NaN; rankErr = NaN;
        return;
      end
      pointsEvaluated = sum(obj.pop.origEvaled);
      pointsToReevaluate = obj.archive.X(end-pointsEvaluated+1:end, :);
      yOriginal = obj.archive.y(end-pointsEvaluated+1:end);
      yReevaled = obj.model.predict(pointsToReevaluate);
      
      % for RMSE and Kendall, get all original-evaled y's except presampled
      %isOriginalAfterPresample = obj.pop.origEvaled & obj.pop.phase ~= 0;
      %yOriginal = obj.pop.y(isOriginalAfterPresample);
      %yReevaled = yModel1(isOriginalAfterPresample);
      rmse = sqrt(sum((yReevaled - yOriginal).^2))/length(yOriginal);
      kendall = corr(yReevaled, yOriginal, 'type', 'Kendall');

      % for RDE, get all (1) the population w/o presampled, and compare it
      % with (2) these points predicted by the first model
      %yModel1 = obj.model.modelPredict(obj.pop.X');
      %yModelMixed = yModel1;
      %yModelMixed(isOriginalAfterPresample) = yOriginal;
      %rankErr = errRankMu(yModel1(obj.pop.phase ~= 0), yModelMixed(obj.pop.phase ~= 0), ...
      %    obj.cmaesState.mu);
      rankErr = NaN;

      % Debug:
      % fprintf('  model: %d pSmpls, reeval %d pts, RMSE= %.2e, Kendl= %.2f, rankErr= %.3f\n', ...
      %     obj.nPresampledPoints, length(yReeval), rmse, kendall, rankErr);
    end


    function [rmse, kCorr, normKendall, age, nData] = oldModelStatistics(obj)
      % return statistics of a several-generations-old model measured
      % on new data from following generations (generations after the
      % model was trained)
      rmse = NaN; kCorr = NaN; normKendall = NaN; age = NaN; nData = 0;

      [oldModel, age] = obj.getOldModel(obj.oldModelAgeForStatistics);
      if (~isempty(oldModel))
        countiter = obj.cmaesState.countiter;
        oldModelGeneration = oldModel.trainGeneration;
        testGenerations = [(oldModel.trainGeneration+1):countiter];
        [xOrig, yOrig] = obj.archive.getDataFromGenerations(testGenerations);
        nData = length(yOrig);
        if (nData > 0)
          yPredict = oldModel.predict(xOrig);
          rmse = sqrt(sum((yPredict - yOrig).^2))/length(yPredict);
          if (nData >= 2)
            kCorr = corr(yPredict, yOrig, 'type', 'Kendall');
            normKendall = (-kCorr + 1) / 2;
          end
        end
      end
    end

    function rankErr = retrainStatistics(obj, yModel1)
    % get ranking error between the first and the second model
    % (if the second model is trained)
      if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
        yModel2 = obj.retrainedModel.predict(obj.pop.x');
        rankErr = errRankMu(yModel1(obj.pop.phase ~= 0), yModel2(obj.pop.phase ~= 0), ...
            obj.cmaesState.mu);
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

      rmse = NaN; kendall = NaN; errRank = NaN; lastModel = [];
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

      % do we have access to the original BBOB fitness?
      if (~ isfield(obj.surrogateOpts.modelOpts, 'bbob_func'))
        return;
      end

      N_VALIDATION_CYCLES = 1;
      kendall = NaN(1, N_VALIDATION_CYCLES);
      rmse = NaN(1, N_VALIDATION_CYCLES);
      errRank = NaN(1, N_VALIDATION_CYCLES);

      for i = 1:N_VALIDATION_CYCLES

      [~, xValidTest, ~] = ...
          sampleCmaesNoFitness(obj.cmaesState.sigma, obj.cmaesState.lambda, obj.cmaesState, sampleOpts);
      preciseModel = ModelFactory.createModel('bbob', obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');
      yTest = preciseModel.predict(xValidTest');
      yPredict = lastModel.predict(xValidTest');

      kendall(i) = corr(yPredict, yTest, 'type', 'Kendall');
      rmse(i) = sqrt(sum((yPredict - yTest).^2))/length(yPredict);

      errRank(i) = errRankMu(yTest, yPredict, obj.cmaesState.mu);

      end

      kendall = nanmean(kendall);
      rmse = nanmean(rmse);
      errRank = nanmean(errRank);

      % Debug:
      % fprintf('  test RMSE= %.2e, Kendall= %.3f, rankErr= %.3f %s\n', ...
      %     rmse, kendall, errRank, decorateKendall(kendall));
    end


    function [obj, yBestOrig, xBest, xBestValid, zBest] = preselection(obj, nBestPoints, nPoints, sampleOpts, varargin)
    % Preselection: orig-evaluate the best predicted point(s)
    % out of obj.preselectionPopRatio*lambda sampled points

      % determine the number of presampled 'best' points
      if (nPoints <= 1 || length(nBestPoints) == 1)
        % if there _is not_ enough points to orig-evalute, use the preselection with
        % different probability...
        obj.usedBestPoints = getProbNumber(nBestPoints(1));
      else
        % ... that if there _is_ enough points
        obj.usedBestPoints = getProbNumber(nBestPoints(2));
      end

      if (nPoints >= 1 && obj.usedBestPoints >= 1)
        % The preselection should really happen.
        % select the points:
        obj.usedBestPoints = min(obj.usedBestPoints, nPoints);
        [~, xBest, xBestValid, zBest] = preselect(obj.usedBestPoints, obj.cmaesState, ...
            obj.model, sampleOpts, obj.preselectionPopRatio);
        % use the original fitness function for them
        [yBestOrig,  xBest, xBestValid, zBest, obj.counteval] = ...
            sampleCmaesOnlyFitness(xBest, xBestValid, zBest, ...
            obj.cmaesState.sigma, obj.usedBestPoints, obj.counteval, obj.cmaesState, sampleOpts, ...
            'Archive', obj.archive, varargin{:});
        % and save them into the Population
        phase = 1;        % this is the first orig-evaluated point
        % remove this number of non-evaluated points
        %
        % TODO: remove not the first points (which are effectively random points),
        %       but be a bit more clever! e.g. the worst predicted points...?
        %
        obj.pop = obj.pop.removeNotOrigEvaluated(obj.usedBestPoints);
        obj.pop = obj.pop.addPoints(xBestValid, yBestOrig, xBest, zBest, obj.usedBestPoints, phase);
        obj.archive.save(xBestValid', yBestOrig', obj.cmaesState.countiter);
      else
        % no preselection
        yBestOrig = []; xBest = []; xBestValid = []; zBest = [];
        obj.usedBestPoints = 0;
      end
    end

    function newPop = updateNotOrigValues(obj, model, pop, phase)
      % updates f-values in the 'pop' with the predictions of the 'model'
      notOrigEvaledX = pop.getNotOrigEvaledX();
      yModel = model.predict(notOrigEvaledX');
      newPop = pop.updateYValue(notOrigEvaledX, yModel', 0, phase, ~pop.origEvaled);
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

function probNum = getProbNumber(exactNumber)
% Calculates randomized non-negative integer as follows:
%   probNum = floor(exactNumber) + eps,
% where eps is 0 or 1. Probability that eps is 1 is equal to the remainder:
%   P[eps = 1] = exactNumber - floor(exactNumber).
  fracN = exactNumber - floor(exactNumber);
  plus = 0;
  if (fracN > 0)
    plus = (rand() < fracN);
  end
  probNum = floor(exactNumber) + plus;
end
