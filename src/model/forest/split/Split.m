classdef Split
% Split creates a split function used in decision trees
  
  properties (Constant)
    splitCandidate = struct( ...
      'splitter', @(X) zeros(size(X, 1), 1), ...
      'gain', -inf ...
      );
  end
  
  properties %(Access = protected)
    transformationOptions % transform options
    transformation % transformation
    X % input data
    y % output data
    allEqual % whether all y values are equal
    soft % use soft split
    lambda % lambda steepness in soft logit function
  end
  
  methods
    function obj = Split(options)
      obj.transformationOptions = defopts(options, 'transformationOptions', struct);
      obj.soft = defopts(options, 'soft', false);
      obj.lambda = defopts(options, 'lambda', 1);
    end
    
    function obj = reset(obj, X, y)
    % sets new transformed input
      [obj.X, obj.y, obj.transformation] = ...
        transform(X, y, obj.transformationOptions);
      n = size(X, 1);
      varY = var(y);
      obj.allEqual = any(isnan(varY)) || all(varY < eps(max(varY)) * n);
    end
    
    function best = get(obj, splitGain)
    % returns the split with max splitGain
      best = obj.splitCandidate;
      if obj.allEqual
        return
      end
      candidate = obj.splitCandidate;
      candidate.gain = splitGain(candidate.splitter);
      best = candidate;
    end
  end
  
  methods (Access = protected)
    function f = createSplitter(obj, splitter)
      trans = obj.transformation;
      if obj.soft
        % soft splitter using logit
        % normalize first
        r = splitter(obj.X);
        idx = r <= 0;
        mm = [mean(-r(idx)), mean(r(~idx))];
        splitter = @(X) Split.normalize(splitter(X), mm);

        lambda = obj.lambda;
        f = @(X) 1 ./ (1 + exp(-lambda * splitter(transformApply(X, trans))));
      else
        % hard splitter
        f = @(X) splitter(transformApply(X, trans)) <= 0;
      end
    end
    
    function f = createModelSplitter(obj, model)
      trans = obj.transformation;
      if obj.soft
        % soft splitter
        f = @(X) Split.modelProbability(model, transformApply(X, trans));
      else
        % hard splitter
        f = @(X) model.predict(transformApply(X, trans)) == 1;
      end
    end
  end
  
  methods (Access = private, Static)
    function p = modelProbability(model, X)
      [~, p] = model.predict(X);
    end
    
    function r = normalize(r, mm)
      idx = r <= 0;
      if mm(1) > 0
        r(idx) = r(idx) / mm(1);
      end
      if mm(2) > 0
        r(~idx) = r(~idx) / mm(2);
      end
    end
  end
end