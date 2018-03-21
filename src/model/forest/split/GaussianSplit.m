classdef GaussianSplit < RandomSplit
% GaussianSplit uses the EM algorithm for Gaussian mixtures to find two 
% clusters in the data and locally transforms the regression problem into 
% a classification problem based on closeness to these clusters.
% Splitting hyperplane provided by LDA or QDA.
%
% Dobra, A. and Gehrke J.: SECRET: a scalable linear regression tree
% algorithm
  
  properties %(Access = protected)
    split_discrimType % degree for discriminant analysis ('linear', 'quadratic')
    split_includeInput % whether to include input space
  end
  
  methods
    function obj = GaussianSplit(options)
      obj = obj@RandomSplit(options);
      obj.split_discrimType = defopts(options, 'split_discrimType', {'linear', 'quadratic'});
      obj.split_includeInput = defopts(options, 'split_includeInput', true);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      [n, d] = size(obj.split_X);
      % all equal or too few data to model Gaussian mixture 
      if obj.split_allEqual || n < d + 2
        return
      end
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
      % repeat loop
      for iRepeats = 1:nRepeats
        % fit 2 clusters in input-output space
        if obj.split_includeInput
          % select random features
          % the amount of features increases with iRepeats
          nFeatures = ceil(d * iRepeats / nRepeats);
          features = datasample(1:d, nFeatures, 'Replace', false);
          Z = [ZX(:, features) Zy];
        else
          Z = Zy;
        end
        warning('off', 'stats:gmdistribution:FailedToConverge');
        model = fitgmdist(Z, 2, 'RegularizationValue', 0.001);
        warning('on', 'stats:gmdistribution:FailedToConverge');
        c = model.cluster(Z);
        
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
