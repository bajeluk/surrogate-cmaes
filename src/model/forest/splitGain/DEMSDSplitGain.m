classdef DEMSDSplitGain < SplitGain
% UMVUESplitGain evaluates split functions used in decision trees using
% uniformly minimum-variance unbiased estimator.
%
% Ahmed (1989): Entropy expressions and their estimators for multivariate
% distributions.
%
% Warning: This split gain function presumes normal distribution. This may
% not hold in case of non-linear regression in leaf, e.g. quadratic
% regression has Chi-squared distribution.

  methods
    function obj = DEMSDSplitGain(options)
      obj = obj@SplitGain(options);
      obj.splitGain_computeSd2 = true;
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, d] = size(data.y);
      % TODO: There should not be an absolute value. However, the
      % determinant is sometimes negative.
      value = d*log(exp(1)*pi) ...
        + log(abs(det(data.y*data.y')) + realmin);
      for i = 1:d
        value = value - digamma( (n + 1 - i) / 2);
      end
      value = 0.5 * value;
    end
  end
end