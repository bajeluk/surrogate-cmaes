classdef GradientSplitGain < SplitGain
% GradientSplitGain evaluates split functions used in decision trees using
% 1st and 2nd gradients of the loss function

  properties
    regularization
  end

  methods
    function obj = GradientSplitGain(options)
      obj = obj@SplitGain(options);
      obj.regularization = defopts(options, 'regularization', 0);
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      % first derivatives
      G = sum(data.y(:, 1));
      % second derivatives
      H = sum(data.y(:, 2));
      value = -0.5 * G*G / (H + obj.regularization);
      % normalize
      value = value / (obj.yRange^2);
    end
  end
  
end