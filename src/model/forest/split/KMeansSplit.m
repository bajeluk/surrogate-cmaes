classdef KMeansSplit < RandomSplit
  
  properties %(Access = protected)
    split_discrimType % degree for discriminant analysis ('linear', 'quadratic')
    split_includeInput % whether to include input space
    split_kmeans_metric % k-means metric (according to matlab kmeans)
    split_kmeans_k      % k-means k parameter
  end
  
  methods
    function obj = KMeansSplit(options)
      obj = obj@RandomSplit(options);
      obj.split_discrimType = defopts(options, 'split_discrimType', {'linear', 'quadratic'});
      obj.split_includeInput = defopts(options, 'split_includeInput', true);
      obj.split_kmeans_metric = defopts(options, 'split_kmeans_metric', 'sqeuclidean');
      obj.split_kmeans_k      = defopts(options, 'split_kmeans_k', 2);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.split_allEqual
        return
      end
      [~, d] = size(obj.split_X);
            % check types of discriminant analysis
      if iscell(obj.split_discrimType)
        discrimTypes = obj.split_discrimType;
      else
        discrimTypes = {obj.split_discrimType};
      end
      % get number of discriminant analysis types, number of repetitions,
      % and number of hyperplanes in last repetition; higher number of
      % repetitions than dimension is useless
      nDiscrTypes = numel(discrimTypes);
      [nRepeats, maxHypRem] = obj.getRepeats(nDiscrTypes, d);
      
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
        c = kmeans(Z, obj.split_kmeans_k, 'Distance', obj.split_kmeans_metric);
        
        % in last repetition check the number of remaining hyperplanes
        if iRepeats == nRepeats
          discrimTypes = discrimTypes(randperm(nDiscrTypes, maxHypRem));
        end
        % discriminant analysis of two clusters
        candidate = obj.getDiscrAnal(splitGain, c, best, discrimTypes);
        if candidate.gain > best.gain
          best = candidate;
        end
        
      end
    end
  end
end