classdef KMeansSplit < RandomSplit
  
  properties %(Access = protected)
    split_discrimType % degree for discriminant analysis ('linear', 'quadratic')
    split_includeInput % whether to include input space
    split_kmeans_metric % k-means metric (according to matlab kmeans)
  end
  
  methods
    function obj = KMeansSplit(options)
      obj = obj@RandomSplit(options);
      obj.split_discrimType = defopts(options, 'split_discrimType', {'linear', 'quadratic'});
      obj.split_includeInput = defopts(options, 'split_includeInput', true);
      obj.split_kmeans_metric = defopts(options, 'split_kmeans_metric', 'sqeuclidean');
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [~, d] = size(obj.split_X);
      % clusters are fit in scaled input-output space
      ZX = zscore(obj.split_X);
      Zy = zscore(obj.split_y);
      for iRepeats = 1:obj.split_nRepeats
        % fit 2 clusters in input-output space
        if obj.split_includeInput
          % select random features
          % the amount of features increases with iRepeats
          nFeatures = ceil(d * iRepeats / obj.split_nRepeats);
          features = datasample(1:d, nFeatures, 'Replace', false);
          Z = [ZX(:, features) Zy];
        else
          Z = Zy;
        end
        c = kmeans(Z, 2, 'Distance', obj.split_kmeans_metric);
        
        % discriminant analysis of two clusters
        if iscell(obj.split_discrimType)
          discrimTypes = obj.split_discrimType;
        else
          discrimTypes = {obj.split_discrimType};
        end
        best = obj.getDiscrAnal(splitGain, c, best, discrimTypes);
        
      end
    end
  end
end