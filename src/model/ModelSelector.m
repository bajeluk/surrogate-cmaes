classdef ModelSelector < Model
  %MODELSELECTOR An abstract composite class aggregating a finite set of models.
  %   Each time training is called, all models are evaluated.
  %   The prediction is performed with the best model in the set.

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

    % model selector-specific properties
    options
    modelOptions
    sharedModelOptions
    factory
    xMean
    models
    modelNames
    modelTypes
    modelIsTrained
    nTrainErrors
    nModels
    bestIdx
    trainTrial
    trainLikelihood
  end

  methods (Access = protected)
    function obj = ModelSelector(options, xMean)
      % comply with Model interface
      obj.xMean = xMean;
      obj.options = options;
      obj.transformCoordinates = defopts(options, 'transformCoordinates', true);
      obj.dim = size(xMean, 2);
      obj.predictionType = defopts(options, 'predictionType', 'fvalues');

      obj.modelOptions = defopts(options, 'modelOptions', []);
      if isempty(obj.modelOptions) || ~isstruct(obj.modelOptions)
        error('ModelSelector: Option ''modelOptions'' must be a struct array with fields ''name'', ''type'' and ''params''');
      end
      obj.nModels = length(obj.modelOptions);
      obj.sharedModelOptions = defopts(options, 'sharedModelOptions', struct());
      if (~isstruct(obj.sharedModelOptions))
        error('ModelSelector: Option ''sharedModelOptions'' must be a struct');
      end

      factory = defopts(options, 'factory', 'ModelFactory');
      if isa(factory, 'function_handle') || ischar(factory)
        obj.factory = feval(factory);
      else
        obj.factory = factory;
      end

      obj.models = cell(1, obj.nModels);
      obj.modelNames = cell(1, obj.nModels);
      obj.modelTypes = cell(1, obj.nModels);
      obj.modelIsTrained = false(1, obj.nModels);
      obj.nTrainErrors = inf(1, obj.nModels);
      obj.trainTrial = 1;

      % merge shared options into specialized options
      % assumes disjunct sets of fields
      for i = 1:obj.nModels
        modelOpts = obj.modelOptions(i).params;
        S = [fieldnames(modelOpts)' fieldnames(obj.sharedModelOptions)'; ...
             struct2cell(modelOpts)' struct2cell(obj.sharedModelOptions)'];
        try
          obj.modelOptions(i).params = struct(S{:});
        catch err
          error('ModelSelector: Could not merge shared options with model-specific options');
        end

        params = obj.modelOptions(i).params;
        name = obj.modelOptions(i).name;
        mdlType = obj.modelOptions(i).type;

        mdl = obj.factory.createModel(mdlType, params, obj.xMean);
        obj.models{i} = mdl;
        obj.modelNames{i} = name;
        obj.modelTypes{i} = mdlType;
      end

      obj.bestIdx = [];
    end

  end

  methods (Abstract)
    [obj, mdlIdx, val] = modelSelect(obj, generation)
    [obj, mdlIdx] = modelSort(obj)
  end

  methods (Access = public)
    function obj = trainModel(obj, X, y, xMean, generation)
      if (~isempty(X) && ~isempty(y))
        obj.dataset.X = X;
        obj.dataset.y = y;
      end

      for i = 1:obj.nModels
        try
          fprintf('Training model %s\n', obj.modelNames{i});
          mdl = obj.models{i}.trainModel(X, y, xMean, generation);
          obj.models{i} = mdl;
          obj.modelIsTrained(obj.trainTrial, i) = obj.models{i}.isTrained();
          obj.nTrainErrors(obj.trainTrial, i) = obj.models{i}.nErrors;
        catch err
          warning('Could not train model %d with error:\n%s\n', i, getReport(err));
          obj.modelIsTrained(obj.trainTrial, i) = false;
        end
      end

      [obj, mdlIdx, val] = obj.modelSelect(obj.trainTrial);

      if any(obj.modelIsTrained(obj.trainTrial, :))
        if isinf(val) || isnan(val)
          mdlIdx = find(obj.modelIsTrained(obj.trainTrial, :), 1);
          warning(['ModelSelector: Invalid IC value for the selected model.\n' ...
            'Falling back to the first model having been trained (%d / %s)'], mdlIdx, obj.modelNames{mdlIdx});
        end
        obj.bestIdx(end+1) = mdlIdx;
        obj.trainGeneration = obj.models{mdlIdx}.trainGeneration;
        obj.trainLikelihood = obj.models{mdlIdx}.trainLikelihood;
      else
        obj.trainGeneration = -1;
      end

      obj.trainTrial = obj.trainTrial + 1;
    end

    function [y, sd2] = modelPredict(obj, X)
      if isempty(obj.bestIdx) || obj.bestIdx(end) == 0
        error('ModelSelector: Best model has not been determined. Was trainModel called?\n');
      end

      [y, sd2] = obj.models{obj.bestIdx(end)}.modelPredict(X);

      % if the best model fails, try the rest in order of their fitness
      if isempty(y) || isempty(sd2)
        idx = obj.modelSort();
        for i = idx(2:end)
          [y, sd2] = obj.models{i}.modelPredict(X);
          if ~isempty(y) && ~isempty(sd2)
            fprintf('ModelSelector: Prediction with model (%d / %s)\n', ...
              i, obj.modelNames{i});
            break;
          end
        end
      end
    end

    function nData = getNTrainData(obj)
      nData = 3 * obj.dim;
    end

  end

end
