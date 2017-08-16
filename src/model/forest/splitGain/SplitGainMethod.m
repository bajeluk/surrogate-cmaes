classdef SplitGainMethod
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
  enumeration
    Data, PolynomialModelModel, ProbabiliticModel, 
  end
end