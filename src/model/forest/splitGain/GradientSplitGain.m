classdef GradientSplitGain < SplitGain
% GradientSplitGain evaluates split functions used in decision trees using
% 1st and 2nd gradients of the loss function

  properties
    splitGain_regularization % regularization
  end

  methods
    function obj = GradientSplitGain(options)
      obj = obj@SplitGain(options);
      obj.splitGain_regularization = defopts(options, 'splitGain_regularization', 0);
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      GH = sum(data.y);
      % first derivatives
      G = GH(1);
      % second derivatives
      H = GH(2);
      value = -0.5 * G*G / (H + obj.splitGain_regularization);
    end
  end
  
end