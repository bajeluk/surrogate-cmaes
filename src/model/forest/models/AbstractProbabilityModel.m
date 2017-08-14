classdef (Abstract) AbstractProbabilityModel < AbstractModel
  
  methods (Access = protected)
    function obj = AbstractProbabilityModel(modelOptions, xMean)
      obj = obj@AbstractModel(modelOptions, xMean);
    end
  end
  
  methods (Abstract)
    [y, sd2, ci, p] = modelPredict(obj, X)
    % predicts the function values in new points X
    % returns empty @y on error
  end
end

