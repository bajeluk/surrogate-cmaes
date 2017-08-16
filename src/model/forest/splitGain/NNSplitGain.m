classdef NNSplitGain < SplitGain
% NNSplitGain evaluates split functions used in decision trees using
% 1 nearest neighbor

  methods
    function obj = NNSplitGain()
    % Creates a new splitter based on 1-NN algorithm
      obj = obj@SplitGain();
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [y,idx] = unique(data.y);
      n = numel(y);
      if n == 1
        value = -realmax / 2;
        return;
      end
      n = numel(y);
      dist = abs(y(1:n-1) - y(2:n));
      nearest = min([dist; inf], [inf; dist]);
      counts = zeros(size(y));
      for i = 1:numel(y)
        counts(idx(i)) = counts(idx(i)) + 1;
      end
      eulergamma = 0.5772;
      value = sum(counts .* log(nearest)) / n + log(n-1) + eulergamma + log(pi/2);
    end
  end
end