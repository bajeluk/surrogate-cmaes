classdef ICModelSelector < ModelSelector
  %
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
    nModels
    bestIdx

    % ICModelSelector properties
    ic % information criterion to use
    modelsIC
  end
  
  methods (Access = protected)
    function calcICs(obj, generation)
      obj.modelsIC.aic(generation, :) = inf;
      obj.modelsIC.bic(generation, :) = inf;

      for mdlIdx = 1:obj.nModels
        if ~obj.isTrained(generation, mdlIdx)
          continue;
        end

        mdl = obj.models{mdlIdx};

        lik = -mdl.getNegLogML();
        k = mdl.getNParams();
        n = mdl.getNData();

        assert(n > 0);

        aic = -2 * lik + 2 * k;
        bic = -2 * lik + log(n) * k;

        obj.modelsIC.aic(generation, mdlIdx) = aic;
        obj.modelsIC.bic(generation, mdlIdx) = bic;
      end
    end
  end

  methods
    function obj = ICModelSelector(modelOptions, xMean)
      obj = obj@ModelSelector(modelOptions, xMean);

      obj.ic = defopts(modelOptions, 'ic', 'bic');

      if ~ismember(obj.ic, {'aic', 'bic'})
        error('Unknown information criterion: ''%s''', obj.ic);
      end

      obj.modelsIC = struct( ...
        'aic', inf(1, obj.nModels), ...
        'bic', inf(1, obj.nModels) ...
      );
    end

    function [mdlIdx, ic] = modelSelect(obj, generation)
      obj.calcICs(generation);
      ics = obj.modelsIC.(obj.ic);
      [ic, mdlIdx] = min(ics(generation, :));
    end
  end

end

