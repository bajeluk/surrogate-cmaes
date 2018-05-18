classdef BayesianICModel < ICModel
  %BAYESIANSELECTABLE Summary of this class goes here
  %   Detailed explanation goes here
  
  methods (Abstract)
    l = getNegLogEst(obj)    % negative log likelihood at a Bayes estimate
    p = getLogPredDens(obj)  % pointwise pointwise predictive density
  end

end

