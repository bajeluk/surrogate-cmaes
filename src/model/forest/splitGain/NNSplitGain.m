classdef NNSplitGain < SplitGain
% NNSplitGain evaluates split functions used in decision trees using
% 1 nearest neighbor

  properties (Access = private)
    sorted = struct('y', [], 'idx', []) % sorted y values
  end
  
  methods
    function obj = NNSplitGain()
    % Creates a new splitter based on 1-NN algorithm
      obj = obj@SplitGain();
    end
  end
  
  methods
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      obj = reset@SplitGain(obj, X, y);
      [obj.sorted.y, obj.sorted.idx] = sort(y);
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      % y = sort(data.y);
      y = nan(size(data.y));
      iNext = 1;
      % loop y values in ascending order
      for i = 1:numel(obj.y)
        j = obj.sortedIdx(i);
        % check if current y is in data.y
        if data.idx(j)
          y(iNext) = obj.y(j);
          iNext = iNext + 1;
        end
      end
      n = numel(y);
      dist = abs(y(1:n-1) - y(2:n));
      nearest = min([dist inf], [inf dist]);
      eulergamma = 0.5772;
      value = sum(log(nearest)) / n + log(n-1) + eulergamma + log(pi/2);
    end
  end
end