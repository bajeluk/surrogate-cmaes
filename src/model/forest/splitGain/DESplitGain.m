classdef DESplitGain < SplitGain
% DESplitGain evaluates split functions used in decision 
% trees using estimate of differential entropy.
% 
% Criminisi(2011): Decision Forests for Classification, Regression, Density
% Estimation, Manifold Learning and Semi-Supervised Learning
%
% Warning: This split gain function presumes normal distribution. This may
% not hold in case of non-linear regression in leaf, e.g. quadratic
% regression has Chi-squared distribution.

  methods
    function obj = DESplitGain(options)
      obj = obj@SplitGain(options);
      obj.splitGain_computeSd2 = true;
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [n, ~] = size(data.y);
      % TODO: data.sd2 is probably not the same as \sigma_y(x) in Criminisi
      value = 1 + log(2*pi) + sum(log(data.sd2 + realmin)) / (2*n);
    end
  end
end