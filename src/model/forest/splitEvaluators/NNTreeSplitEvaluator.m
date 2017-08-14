classdef NNTreeSplitEvaluator < TreeSplitEvaluator
  
  properties (Constant, Access = private)
    nSamples = 256;
  end
  
  properties (Access = private)
    metric % distance metric
    D % distances
  end
  
  methods
    function obj = NNTreeSplitEvaluator(metric, sample)
      obj = obj@TreeSplitEvaluator();
      obj.metric = metric;
      obj.sampe = sample;
    end
  end
  
  methods (Access = protected)
    function resetData(obj)
      if ~obj.sample
        D = squareform(pdist(y, obj.metric));
        obj.D = D + inf * eye(size(D));
      end
    end
    
    function value = getValue(obj, data)
      if obj.sample
        n = min(size(data.y, 1), NNTreeSplitEvaluator.nSamples);
        y = datasample(data.y, n, 'Replace', false);
        D = squareform(pdist(y, obj.metric));
        D = D + inf * eye(size(D));
      else
        n = size(data.y, 1);
        D = obj.D(data.idx, data.idx);
      end
      nearestD = min(D)';
      eulergamma = 0.5772;
      value = sum(log(nearestD)) / n + log(n-1) + eulergamma + log(pi/2);
    end
  end
end