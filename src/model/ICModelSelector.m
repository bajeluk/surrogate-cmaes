classdef ICModelSelector < ModelSelector
  %
  properties
    % ICModelSelector properties
    ic % information criterion to use
    modelsIC
  end
  
  methods (Access = protected)
    function obj = calcICs(obj, generation)
      obj.modelsIC.aic(generation, :) = inf(1, obj.nModels);
      obj.modelsIC.bic(generation, :) = inf(1, obj.nModels);

      for mdlIdx = 1:obj.nModels
        if ~obj.modelIsTrained(generation, mdlIdx)
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

      obj.modelsIC = struct( ...
        'aic', inf(1, obj.nModels), ...
        'bic', inf(1, obj.nModels) ...
      );
    end

    function [obj, mdlIdx, ic] = modelSelect(obj, generation)
      obj = obj.calcICs(generation);
      ics = obj.modelsIC.(obj.ic);

      [ic, mdlIdx] = min(ics(generation, :));
    end
  end

end

