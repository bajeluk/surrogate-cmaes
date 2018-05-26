classdef ICModel
  %SELECTABLE A model that supports calculation of information criteria.
  
  properties
  end
  
  methods (Abstract)
    k = getNParams(obj)     % number of parameters in the model
    n = getNData(obj)       % number of data points in the model
    l = getNegLogML(obj)    % estimation of negative log likelihood mode
    nloo = getNegLooPredDens(obj)    % negative leave-one-out predictive density
  end

end

