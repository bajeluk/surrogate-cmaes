classdef BayesianICModel < ICModel
  %BAYESIANSELECTABLE Summary of this class goes here
  %   Detailed explanation goes here
  
  methods (Abstract)
    l = getNegLogEst(obj)              % negative log likelihood at a Bayes estimate
    p = getLogPredDens(obj, varargin)  % pointwise pointwise predictive density
    s = getChains(obj)                 % MCMC chains
    lik = getNegLogLPost(obj)          % likelihoods of a sample
  end

end

