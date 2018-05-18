classdef GPICModelSelector < BayesianICModelSelector
  methods (Access = public)
    function obj = GPICModelSelector(modelOptions, xMean, options)
      obj = obj@BayesianICModelSelector(modelOptions, xMean, options);
    end
  end

  methods (Access = protected)
    function createModel(obj, mdlIdx)
      params = obj.modelOptions(mdlIdx).params;
      name = obj.modelOptions(mdlIdx).name;

      mdl = GpModel(params, obj.xMean);
      obj.models{mdlIdx} = mdl;
      obj.modelNames{mdlIdx} = name;
    end
  end
end

