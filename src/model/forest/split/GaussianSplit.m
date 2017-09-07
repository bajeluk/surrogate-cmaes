classdef GaussianSplit < RandomSplit
  
  properties %(Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    includeInput % whether to include input space
  end
  
  methods
    function obj = GaussianSplit(options)
      obj = obj@RandomSplit(options);
      obj.discrimType = defopts(options, 'discrimType', {'linear', 'quadratic'});
      obj.includeInput = defopts(options, 'includeInput', true);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      trans = obj.transformation;
      [n, d] = size(obj.X);
      % clusters are fit in scaled input-output space
      ZX = zscore(obj.X);
      Zy = zscore(obj.y);
      for iRepeats = 1:obj.nRepeats
        candidate = obj.splitCandidate;
        % fit 2 clusters in input-output space
        if obj.includeInput
          % select random features
          % the amount of features increases with iRepeats
          nFeatures = ceil(d * iRepeats / obj.nRepeats);
          features = datasample(1:d, nFeatures, 'Replace', false);
          Z = [ZX(:, features) Zy];
        else
          Z = Zy;
        end
        warning('off', 'stats:gmdistribution:FailedToConverge');
        model = fitgmdist(Z, 2, 'RegularizationValue', 0.001);
        warning('on', 'stats:gmdistribution:FailedToConverge');
        c = model.cluster(Z);
        
        % discriminant analysis of two clusters
        if iscell(obj.discrimType)
          discrimTypes = obj.discrimType;
        else
          discrimTypes = {obj.discrimType};
        end
        for i = 1:numel(discrimTypes)
          discrimType = discrimTypes{i};
          try
            model = fitcdiscr(obj.X, c, 'DiscrimType', discrimType);
          catch
            % singular covariance matrix
            pseudoDiscrimType = strcat('pseudo', discrimType);
            model = fitcdiscr(obj.X, c, 'DiscrimType', pseudoDiscrimType);
          end
          candidate.splitter = @(X)...
            model.predict(transformApply(X, trans)) == 1;
          candidate.gain = splitGain.get(candidate.splitter);
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
    end
  end
  
end