classdef ModelSelector < Model
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
    nModels
    bestIdx
  end

  methods (Access = protected)
    function obj = ModelSelector(options, xMean)
      obj.xMean = xMean;
      obj.selectorOptions = options;
      obj.modelOptions = defopts(options, 'modelOptions', []);
      if ~isempty(obj.models) || ~isstruct(obj.models)
        error('ModelSelector: Option models must be a struct array with fields ''name'' and ''params''');
      end
      obj.nModels = length(obj.modelOptions);
      obj.sharedModelOpts = defopts(options, 'sharedModelOptions', struct());
      if (~isstruct(obj.sharedModelOpts))
        error('ModelSelector: Option ''sharedModelOpts'' must be a struct');
      end
      
      obj.models = cell(1, obj.nModels);

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
          obj.models{i} = obj.createModel(i);
        catch err
          error('ModelSelector: Could not create model %d, %s', ...
            i, obj.modelOptions(i).name);
        end
      end

      obj.bestIdx = 0;
    end
  end

  methods (Abstract)
    mdl = createModel(obj, mdlIdx)
    mdlIdx = modelSelect(obj)
  end

  methods (Access = public)
    function obj = trainModel(obj, X, y, xMean, generation)
      for i = 1:obj.nModels
        obj.models{i}.trainModel(X, y, xMean, generation);
      end

      obj.bestIdx = obj.modelSelect();
    end

    function [y, sd2] = modelPredict(obj, X)
    end
  end

end
