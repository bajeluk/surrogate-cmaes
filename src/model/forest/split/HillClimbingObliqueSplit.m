classdef HillClimbingObliqueSplit < RandomSplit
% HillClimbingObliqueSplit starts with a random hyperplane (or with some
% good axis parallel hyperplane) and deterministically perturbates
% the hyperplane's direction in each axis to maximize the split gain. Once
% no improvement is possible, performs a number of random jumps as an
% attempt to escape local maxima. If some random jump succeeds,
% deterministic perturbation is performed again.
  
  properties
    X1 % X with intercept
    nQuantize % quantization of tresholds
    nRandomPerturbations % number of random perturbations
  end

  methods
    function obj = HillClimbingObliqueSplit(options)
      obj = obj@RandomSplit(options);
      obj.nQuantize = defopts(options, 'nQuantize', 0);
      obj.nRandomPerturbations = defopts(options, 'nRandomPerturbations', 10);
    end
    
    function obj = reset(obj, X, y)
    % sets new transformed input
      obj = reset@RandomSplit(obj, X, y);
      obj.X1 = generateFeatures(obj.X, 'linear', true, true);
    end

    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      for iRepeats = 1:obj.nRepeats
        if iRepeats == 1
          candidate = obj.getAxisHyperplane(splitGain);
        else
          candidate = obj.getRandomHyperplane(splitGain);
        end
        candidate = obj.hillClimb(splitGain, candidate);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
  end
  
  methods (Access = private)
    function best = getAxisHyperplane(obj, splitGain)
      best = obj.splitCandidate;
      trans = obj.transformation;
      [n, d] = size(obj.X);
      for feature = 1:d
        featureSelector = (1:d == feature)';
        values = obj.X(:, feature)';
        if obj.nQuantize > 0 && numel(values) > obj.nQuantize
          mm = minmax(values);
          tresholds = linspace(mm(1), mm(2), obj.nQuantize);
        else
          tresholds = unique(values);
        end
        for treshold = tresholds
          candidate = obj.splitCandidate;
          candidate.splitter = @(X)...
            transformApply(X, trans) * featureSelector <= treshold;
          candidate.gain = splitGain.get(candidate.splitter);
          candidate.feature = feature;
          candidate.treshold = treshold;
          if candidate.gain > best.gain
            best = candidate;
          end
        end
      end
      
      best.H = [zeros(1, d), -best.treshold]';
      best.H(best.feature) = 1;
    end
    
    function candidate = getRandomHyperplane(obj, splitGain)
      [n, d] = size(obj.X);
      H = tan(rand(1, d+1)' * pi - pi/2);
      candidate = obj.getSplit(splitGain, H);
      candidate.H = H;
    end
    
    function best = hillClimb(obj, splitGain, best)
      d = size(obj.X, 2);
      pStag = 0.1;
      pMove = pStag;
      J = obj.nRandomPerturbations;
      while J > 0
        improvement = true;
        while improvement
          improvement = false;
          for feature = 1:d+1
            candidate = obj.deterministicPerturb(splitGain, best, feature);
            if candidate.gain > best.gain
              best = candidate;
              pMove = pStag;
              improvement = true;
            elseif candidate.gain == best.gain
              if rand() < pMove
                best = candidate;
                pMove = pMove - 0.1 * pStag;
              end
            end
          end
        end
        while J > 0
          candidate = obj.randomPerturb(splitGain, best);
          if candidate.gain >= best.gain
            best = candidate;
            break;
          end
          J = J - 1;
        end
      end
    end
    
    function best = deterministicPerturb(obj, splitGain, best, feature)
      V = obj.X1 * best.H;
      U = (best.H(feature) * obj.X1(:, feature) - V) ...
        ./ obj.X1(:, feature);
      H = best.H;
      values = U';
      if obj.nQuantize > 0 && numel(values) > obj.nQuantize
        mm = minmax(values);
        tresholds = linspace(mm(1), mm(2), obj.nQuantize);
      else
        tresholds = unique(values);
      end
      for treshold = tresholds
        if isinf(treshold) || isnan(treshold)
          continue
        end
        H(feature) = treshold;
        candidate = obj.getSplit(splitGain, H);
        if candidate.gain > best.gain
          best = candidate;
        end
      end
    end
    
    function candidate = randomPerturb(obj, splitGain, best)
      r = randn(size(best.H));
      r = max(-pi/2, r);
      r = min(pi/2, r);
      r = tan(r);
      H = best.H + r;
      candidate = obj.getSplit(splitGain, H);
    end
    
    function [candidate] = getSplit(obj, splitGain, H)
      candidate = obj.splitCandidate;
      candidate.splitter = obj.createSplitter(@(X) ...
        generateFeatures(X, 'linear', true, true) * H);
      candidate.gain = splitGain.get(candidate.splitter);
      candidate.H = H;
    end
  end
end