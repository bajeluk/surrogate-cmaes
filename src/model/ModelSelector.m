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
    xMean
    models
    modelNames
    modelTypes
    isTrained
    nModels
    bestIdx
    trainTrial
  end

  methods (Access = protected)
    function obj = ModelSelector(options, xMean)
      obj.xMean = xMean;
      obj.selectorOptions = options;
      obj.modelOptions = defopts(options, 'modelOptions', []);
      if ~isempty(obj.models) || ~isstruct(obj.models)
        error('ModelSelector: Option models must be a struct array with fields ''name'', ''type'' and ''params''');
      end
      obj.nModels = length(obj.modelOptions);
      obj.sharedModelOpts = defopts(options, 'sharedModelOptions', struct());
      if (~isstruct(obj.sharedModelOpts))
        error('ModelSelector: Option ''sharedModelOpts'' must be a struct');
      end

      if ~isfield(options, 'factory') || isempty(options.factory)
        error('ModelSelector: Model factory must be supplied to create models on demand.');
      else
        factory = options.factory
        if isa(factory, 'function_handle') || ischar(factory)
          obj.factory = feval(factory);
        else
          obj.factory = factory;
        end
      end

      obj.models = cell(1, obj.nModels);
      obj.modelNames = cell(1, obj.nModels);
      obj.modelTypes = cell(1, obj.nModels);
      obj.isTrained = false(1, obj.nModels);
      obj.trainTrail = 1;

      % merge shared options into specialized options
      % assumes disjunct sets of fields
      for i = 1:obj.nModels
        modelOpts = obj.modelOptions(i).params;
        S = [fieldnames(modelOpts)' fieldnames(obj.sharedModelOpts)'; ...
             struct2cell(modelOpts)' struct2cell(obj.sharedModelOpts)'];
        try
          obj.modelOptions(i).params = struct(S{:});
        catch err
          error('ModelSelector: Could not merge shared options with model-specific options');
        end

        try
          obj.createModel(i);
        catch err
          error('ModelSelector: Could not create model %d, %s', ...
            i, obj.modelOptions(i).name);
        end
      end

      obj.bestIdx = 0;
    end

    function createModel(obj, mdlIdx)
      params = obj.modelOptions(mdlIdx).params;
      name = obj.modelOptions(mdlIdx).name;
      mdlType = obj.modelOptions(mdlIdx).type;

      mdl = obj.factory.createModel(mdlType, params, obj.xMean);
      obj.models{mdlIdx} = mdl;
      obj.modelNames{mdlIdx} = name;
      obj.modelTypes{mdlIdx} = mdlType;
    end
  end

  methods (Abstract)
    [mdlIdx, val] = modelSelect(obj, generation)
  end

  methods (Access = public)
    function obj = trainModel(obj, X, y, xMean, generation)
      for i = 1:obj.nModels
        try
          obj.models{i} = obj.models{i}.trainModel(X, y, xMean, generation);
          obj.isTrained(obj.trainTrial, i) = obj.models{i}.isTrained();
        catch err
          obj.isTrained(obj.trainTrial, i) = false;
        end
      end

      [mdlIdx, val] = obj.modelSelect(obj.trainTrial);

      if any(obj.isTrained(obj.trainTrial, :))
        if isinf(val) || isnan(val)
          warning('ModelSelector: Invalid IC value for the selected model.');
          obj.trainGeneration = -1;
        else
          obj.bestIdx = mdlIdx;
          obj.trainGeneration = obj.models{mdlIdx}.trainGeneration;
        end
      else
        obj.trainGeneration = -1;
      end

      obj.trainTrial = obj.trainTrial + 1;
    end

    function [y, sd2] = modelPredict(obj, X)
      if obj.bestIdx == 0
        error('ModelSelector: Best model has not been determined. Was trainModel called?');
      end
      [y, sd2] = obj.models{obj.bestIdx}.modelPredict(X);
    end

  end

end
