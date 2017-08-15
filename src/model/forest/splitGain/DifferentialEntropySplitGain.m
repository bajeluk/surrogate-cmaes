classdef DifferentialEntropySplitGain < SplitGain
% DifferentialEntropySplitGain evaluates split functions used in decision 
% trees using differential entropy

  methods
    function obj = DifferentialEntropySplitGain(varargin)
      obj = obj@SplitGain(varargin{:});
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      value = 0.5 * ( ...
        1 - log(2*pi) + sum(log(data.sd2)) ...
        );
    end
  end
end