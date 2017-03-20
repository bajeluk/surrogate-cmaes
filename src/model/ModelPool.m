classdef ModelPool < Model
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

    % ModelPool specific properties

    modelPoolOptions
    archive
    modelsCount           % number of models in the modelPool
    historyLength         % number of older models that are saved
    models                % instances of models, 2D cell array of modelscount*historyLength+1
    isModelTrained        % 2D array, 0 if model at this position in models property is trained, 1 otherwise
    bestModelIndex
    bestModelsHistory     % how many times has been each model chosen as the best one
    bestModelSelection    % which criterium will be used to select which model is the best
                          %(mse/mae/rdeAll/rdeOrig/likelihood)
    choosingCriterium     % values of calculated criterium (mse/mae/...) for each model
    retrainPeriod
    nTrainData            % min of getNTrainData of all created models
    xMean
    minTrainedModelsPercentilForModelChoice % if percentile of oldest models that are trained
                                            % drops below this value, we try newer generations of models
    maxGenerationShiftForModelChoice  % stops trying to find trained generation
                                      % and switches to likelihood after this value of searched generations
  end

  methods (Access = public)
    function obj = ModelPool(modelOptions, xMean)
      obj.modelPoolOptions = modelOptions;
      obj.xMean = xMean;
      obj.modelsCount = length(modelOptions.parameterSets);
      assert(obj.modelsCount ~= 0, 'ModelPool(): No model provided!');

      if (strcmpi(obj.bestModelSelection, 'likelihood'))
        % likelihood selection does not need older models
        obj.historyLength = 0;
      else
        obj.historyLength = defopts(modelOptions, 'historyLength', 4);
        if (obj.historyLength < 2)
          warning('ModelPool: history length needs to be at least 2 in order to choose from at least 1 point, choosing 2 as value.');
          obj.historyLength = 2;
        end
        obj.minTrainedModelsPercentilForModelChoice = defopts(modelOptions, 'minTrainedModelsPercentilForModelChoice', 0.25);
        obj.maxGenerationShiftForModelChoice = defopts(modelOptions, 'maxGenerationShiftForModelChoice', 1);
        if obj.maxGenerationShiftForModelChoice >= obj.historyLength -1
          warning('ModelPool: maxGenerationShiftForModelChoice is too high, choosing %d as value.', obj.historyLength - 2)
          obj.maxGenerationShiftForModelChoice = obj.historyLength - 2;
        end
      end

      obj.retrainPeriod = defopts(modelOptions, 'retrainPeriod', 1);
      obj.models = cell(obj.modelsCount,obj.historyLength+1);
      obj.isModelTrained = false(obj.modelsCount,obj.historyLength+1);
      obj.bestModelsHistory = zeros(1,obj.modelsCount);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      obj.stdY      = 1;
      obj.bestModelSelection = defopts(modelOptions, 'bestModelSelection', 'mse');

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.dimReduction = defopts(modelOptions, 'dimReduction', 1);
      obj.options.normalizeY = defopts(modelOptions, 'normalizeY', true);

      obj.nTrainData = Inf;
      for i=1:obj.modelsCount
        %create the models, calculate needed properties
        modelOptions = obj.modelPoolOptions.parameterSets(i);
        obj.modelPoolOptions.parameterSets(i).calculatedTrainRange = ModelPool.calculateTrainRange(modelOptions.trainRange, obj.dim);
        obj.models{i,1} = obj.createGpModel(i, xMean);
        obj.nTrainData = min(obj.models{i,1}.getNTrainData(),obj.nTrainData);
      end
    end

    function gpModel = createGpModel(obj, modelIndex, xMean)
      newModelOptions = obj.modelPoolOptions.parameterSets(modelIndex);
      newModelOptions.predictionType = obj.predictionType;
      newModelOptions.transformCoordinates = obj.transformCoordinates;
      newModelOptions.dimReduction = obj.dimReduction;
      newModelOptions.options.normalizeY = obj.options.normalizeY;

      gpModel = ModelFactory.createModel('gp', newModelOptions, xMean);

    end

    function nData = getNTrainData(obj)
      nData = obj.nTrainData;
    end

    function trained = isTrained(obj)
      % check whether the model chosen as the best in the newest generation is trained
      if (isempty(obj.isModelTrained(obj.bestModelIndex,1)))
        trained = false;
      else
        trained = obj.isModelTrained(obj.bestModelIndex,1);
      end
    end

    function obj = trainModel(obj, ~, ~, ~, ~)
      % This function is empty because it is not needed, training is
      % done in train().
    end

    function obj = train(obj, X, y, stateVariables, sampleOpts, archive, population)
      obj.archive = archive;
      obj.stateVariables = stateVariables;
      obj.sampleOpts = sampleOpts;
      obj.xMean = stateVariables.xmean';
      generation = stateVariables.countiter;
      if (mod(generation,obj.retrainPeriod)==0)

        trainedModelsCount=0;
        for i=1:obj.modelsCount

          obj.models(i,:) = circshift(obj.models(i,:),[0,1]);
          obj.isModelTrained(i,:) = circshift(obj.isModelTrained(i,:),[0,1]);
          obj.isModelTrained(i,1) = 0;
          obj.models{i,1} = obj.createGpModel(i, obj.xMean);
          obj.models{i,1} = obj.models{i,1}.train(X, y, stateVariables, sampleOpts, obj.archive, population);

          if (obj.models{i,1}.isTrained())
            trainedModelsCount = trainedModelsCount+1;
            obj.isModelTrained(i,1) = 1;
          else
            obj.models{i,1}.trainGeneration = -1;
          end
        end

        if (trainedModelsCount==0)
          warning('ModelPool.trainModel(): trainedModelsCount == 0');
        else
          obj.trainGeneration = generation;

          [obj.bestModelIndex,obj.choosingCriterium] = obj.chooseBestModel(generation, population);
          obj.bestModelsHistory(obj.bestModelIndex) = obj.bestModelsHistory(obj.bestModelIndex)+1;
          obj = obj.copyPropertiesFromBestModel();
        end
      end
    end

    function [y, sd2] = modelPredict(obj, X)
      [y,sd2] = obj.models{obj.bestModelIndex,1}.modelPredict(X);
    end

    function X = getDataset_X(obj)
      X = obj.models{obj.bestModelIndex,1}.getDataset_X();
    end

    function y = getDataset_y(obj)
      y = obj.models{obj.bestModelIndex,1}.getDataset_y();
    end

  end

  methods (Access = public)
    function [bestModelIndex, choosingCriterium] = chooseBestModel(obj, lastGeneration, population)

      if (isempty(lastGeneration))
        lastGeneration = 0;
      end

      ageOfTestedModels = -1;
      for i = obj.historyLength+1 : -1 : obj.historyLength+1 - obj.maxGenerationShiftForModelChoice;
        trainedPercentile = mean(obj.isModelTrained(:,i));

        if (trainedPercentile >= obj.minTrainedModelsPercentilForModelChoice)
          ageOfTestedModels = i;
          break;
        end
      end

      if (ageOfTestedModels == -1 || strcmpi(obj.modelPoolOptions.bestModelSelection,'likelihood'))
        choosingCriterium = obj.getLikelihood();
      else
        switch lower(obj.bestModelSelection)
          case 'rdeorig'
            choosingCriterium = obj.getRdeOrig(ageOfTestedModels, lastGeneration);
          case 'rdeall'
            choosingCriterium = obj.getRdeAll(ageOfTestedModels, population);
          case 'mse'
            choosingCriterium = obj.getMse(ageOfTestedModels, lastGeneration);
          case 'mae'
            choosingCriterium = obj.getMae(ageOfTestedModels, lastGeneration);
          otherwise
            error(['ModelPool.chooseBestModel: ' obj.modelPoolOptions.bestModelSelection ' -- no such option available']);
        end
      end
      % choose the best model from trained ones according to the choosing criterium
      [minValue,bestModelIndex] = min(choosingCriterium(obj.isModelTrained(:,1)));
      if minValue==Inf
        bestModelIndex = 1;
        if (mean(obj.isModelTrained(:,i))>0)
          warning('ModelPool.chooseBestModel: value of minimum is Inf, ageOfTestedModels %d, percentile of trainedModels %d', ...
            ageOfTestedModels, mean(obj.isModelTrained(:,i)));
        end
      end
    end

    function obj = copyPropertiesFromBestModel(obj)
      obj.stdY = obj.models{obj.bestModelIndex,1}.stdY;
      obj.options = obj.models{obj.bestModelIndex,1}.options;
      obj.hyp = obj.models{obj.bestModelIndex,1}.hyp;
      obj.meanFcn = obj.models{obj.bestModelIndex,1}.meanFcn;
      obj.covFcn = obj.models{obj.bestModelIndex,1}.covFcn;
      obj.likFcn = obj.models{obj.bestModelIndex,1}.likFcn;
      obj.infFcn = obj.models{obj.bestModelIndex,1}.infFcn;
      obj.nErrors = obj.models{obj.bestModelIndex,1}.nErrors;
      obj.trainLikelihood = obj.models{obj.bestModelIndex,1}.trainLikelihood;

      obj.shiftY = obj.models{obj.bestModelIndex,1}.shiftY;
      obj.trainSigma = obj.models{obj.bestModelIndex,1}.trainSigma;
      obj.trainBD = obj.models{obj.bestModelIndex,1}.trainBD;
      obj.trainMean = obj.models{obj.bestModelIndex,1}.trainMean;
      obj.shiftMean = obj.models{obj.bestModelIndex,1}.shiftMean;
      obj.reductionMatrix = obj.models{obj.bestModelIndex,1}.reductionMatrix;
      obj.stateVariables = obj.models{obj.bestModelIndex,1}.stateVariables;
      obj.sampleOpts = obj.models{obj.bestModelIndex,1}.sampleOpts;
    end

    function choosingCriterium = getRdeOrig(obj, ageOfTestedModels, lastGeneration)
      choosingCriterium = Inf(obj.modelsCount,1);
      for i=1:obj.modelsCount
        generations = obj.models{i,ageOfTestedModels}.trainGeneration+1:lastGeneration;
        [X,yArchive] = obj.archive.getDataFromGenerations(generations);
        if (size(X,1)~=0)
          if (obj.isModelTrained(i,ageOfTestedModels))
            [yModel, ~] = obj.models{i,ageOfTestedModels}.modelPredict(X);
            if (size(yArchive)==size(yModel))
              choosingCriterium(i) = errRankMu(yModel,yArchive,size(yArchive,1));
            end
          end
        end
      end
    end

    function choosingCriterium = getRdeAll(obj, ageOfTestedModels, population, mu)
      choosingCriterium = Inf(obj.modelsCount,1);
      for modelIndex=1:obj.modelsCount
        partialCriteriumValues = [];
        for modelAge=ageOfTestedModels : -1 : 2
          if obj.isModelTrained(modelIndex,modelAge)
            testedModel = obj.models{modelIndex,modelAge};
            nextModel = obj.models{modelIndex,modelAge-1};

            [~, xSample] = sampleCmaesNoFitness(...
              nextModel.trainSigma, ...
              nextModel.stateVariables.lambda, ...
              nextModel.stateVariables, ...
              nextModel.sampleOpts);
            % get points from archive
            [origPoints_X, origPoints_y] = obj.archive.getDataFromGenerations(testedModel.trainGeneration+1);
            if (modelAge == 2) %we are testing newest model, add points from population if they are available
                origPoints_X = [ origPoints_X ; population.x(:,population.origEvaled)'];
                origPoints_y = [ origPoints_y ; population.y(:,population.origEvaled)'];
            end
            xSample(:,1:size(origPoints_X,1)) = origPoints_X(1:size(origPoints_X,1),:)';
            ySample = testedModel.predict(xSample');
            yWithOrig = ySample;
            % replace predicted values with original values
            yWithOrig(1:size(origPoints_y,1)) = origPoints_y;
            partialCriteriumValues(end+1,1) = errRankMu(ySample, yWithOrig, obj.stateVariables.mu);
          else
            partialCriteriumValues(end+1,1) = NaN;
          end
        end
        x = partialCriteriumValues;
        n = length(x);
        weights = 2.^(-(1:n)) ./ sum(2.^(-(1:n)));
        choosingCriterium(modelIndex) = weights(~isnan(x)) * x(~isnan(x)) ./ sum(weights(~isnan(x)));

      end
    end

    function choosingCriterium = getMse(obj, ageOfTestedModels, lastGeneration)
      choosingCriterium = Inf(obj.modelsCount,1);
      for i=1:obj.modelsCount
        generations=obj.models{i,ageOfTestedModels}.trainGeneration+1:lastGeneration;
        [X,yArchive] = obj.archive.getDataFromGenerations(generations);
        if (size(X,1)~=0)
          if (obj.isModelTrained(i,ageOfTestedModels))
            [yModel, ~] = obj.models{i, ageOfTestedModels}.modelPredict(X);
            if (size(yArchive)==size(yModel))
              choosingCriterium(i) = sum((yModel - yArchive).^2)/size(yArchive);
            end
          end
        end
      end
    end

    function choosingCriterium = getMae(obj, ageOfTestedModels, lastGeneration)
      choosingCriterium = Inf(obj.modelsCount,1);
      for i=1:obj.modelsCount
        generations=obj.models{i,ageOfTestedModels}.trainGeneration+1:lastGeneration;
        [X,yArchive] = obj.archive.getDataFromGenerations(generations);
        if (size(X,1)~=0)
          if (obj.isModelTrained(i,ageOfTestedModels))
            [yModel, ~] = obj.models{i,ageOfTestedModels}.modelPredict(X);
            if (size(yArchive)==size(yModel))
              choosingCriterium(i) = sum(abs(yModel - yArchive))/size(yArchive);
            end
          end
        end
      end
    end

    function choosingCriterium = getLikelihood(obj)
      choosingCriterium = Inf(obj.modelsCount,1);
      for i=1:obj.modelsCount
        choosingCriterium(i) = obj.models{i,1}.trainLikelihood;
      end
    end

  end

  methods (Static)
    function result = calculateTrainRange(percentile, dimension)
      result = chi2inv(percentile,dimension);
    end
  end
end
