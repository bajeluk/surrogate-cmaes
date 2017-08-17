classdef GaussianSplit < Split
  
  properties (Access = protected)
    discrimType % degree for discriminant analysis ('linear', 'quadratic')
    includeInput % whether to include input space
    nRepeats % number or repeats
  end
  
  methods
    function obj = GaussianSplit(...
        transformationOptions, discrimType, includeInput, nRepeats)
      obj = obj@Split(transformationOptions);
      obj.discrimType = discrimType;
      obj.includeInput = includeInput;
      obj.nRepeats = nRepeats;
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      % gaussians are fit on scaled input-output space
      if obj.includeInput
        Z = [X y];
      else
        Z = y;
      end
      Z = zscore(Z);
      for iRepeats = 1:obj.nRepeats
        candidate = obj.splitCandidate;
        % fit 2 gaussian distributions on input-output space
        model = fitgmdist(Z, 2, 'RegularizationValue', 0.001);
        c = model.cluster(X);
        % discriminant analysis of two clusters
        model = fitcdiscr(X, c, 'DiscrimType', obj.discrimType);
        candidate.splitter = @(X)...
          model.predict(transformApply(X, trans)) == 1;
        candidate.gain = splitGain.get(splitter);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
end