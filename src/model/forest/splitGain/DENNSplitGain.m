classdef DENNSplitGain < SplitGain
% NNSplitGain evaluates split functions used in decision trees using
% 1 nearest neighbor

  properties
    k % k-th neighbor
  end

  methods
    function obj = DENNSplitGain(options)
    % Creates a new splitter based on k-NN algorithm
      obj = obj@SplitGain(options);
      obj.k = defopts(options, 'k', 1);
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [yUnique, ~, idxUnique] = unique(data.y);
      n = numel(data.y);
      nUnique = numel(yUnique);
      if nUnique <= obj.k
        value = 0; %-realmax / numel(obj.y) * n;
        return;
      end
      dist = abs(yUnique(1:nUnique-obj.k) - yUnique(obj.k+1:nUnique));
      padding = inf(obj.k, 1);
      nearest = min([dist; padding], [padding; dist]);
      counts = zeros(nUnique, 1);
      for i = 1:n
        counts(idxUnique(i)) = counts(idxUnique(i)) + 1;
      end
      eulergamma = 0.5772;
      value = sum(counts .* log(nearest)) / n + log(n-1) + eulergamma + log(pi/2);
    end
  end
end