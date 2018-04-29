classdef DENNSplitGain < SplitGain
% NNSplitGain evaluates split functions used in decision trees using
% k-th nearest neighbor (in 1D).
%
% Kozachenko (1987): Sample estimate of the entropy of a random vector.
% (english description - Beirlant (2001): Nonparametric entropy estimation:
% An overview)

  properties
    splitGain_k % k-th neighbor
  end

  methods
    function obj = DENNSplitGain(options)
    % Creates a new splitter based on k-NN algorithm
      obj = obj@SplitGain(options);
      obj.splitGain_k = defopts(options, 'splitGain_k', 1);
    end
  end
  
  methods (Access = protected)    
    function value = getValue(obj, data)
    % evaluates data using custom metric
      [yUnique, ~, idxUnique] = unique(data.y);
      n = numel(data.y);
      nUnique = numel(yUnique);
      if nUnique <= obj.splitGain_k
        value = 0; %-realmax / numel(obj.splitGain_y) * n;
        return;
      end
      dist = abs(yUnique(1:nUnique-obj.splitGain_k) - yUnique(obj.splitGain_k+1:nUnique));
      padding = inf(obj.splitGain_k, 1);
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