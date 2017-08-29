classdef (Abstract) SplitGain
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
  
  properties %(Access = protected)
    X % input data
    y % output data
    minSize % min size of one side of the split
    modelFunc % function which creates new model, model = modelFunc(xMean)
    probabilistic % whether to compute also probability and confidence interval
    XP % generated polynomial features for MethodPolynomialModel
    degree % degree of polynomial features for MethodPolynomialModel
    current % struct holding data for current node
    allEqual % whether all y values are equal
  end
  
  methods
    function obj = SplitGain(options)
    % constructor
    % evaluates functions based on constant model
    %   (y, yPred, variancePred)
    % if 'degree' or 'modelFunc' is specified, uses custom model
    %   (y, yPred, variancePred)
    % if also 'probabilistic' is specified
    %   (y, yPred, variancePred, confidenceIntervalPred, probabilityPred)
    % options
    %   'minSize'       - specifies the min size of either side
    %   'degree'        - uses polynomial model of given degree
    %   'modelFunc'     - uses custom model
    %   'probabilistic' - computes probability and confidence interval
      if nargin >= 1
        obj.minSize = defopts(options, 'minSize', 1);
        obj.degree = defopts(options, 'degree', []);
        obj.modelFunc = defopts(options, 'modelFunc', []);
        obj.probabilistic = defopts(options, 'probabilistic', false);
      end
    end
    
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      [n, ~] = size(X);
      obj.X = X;
      obj.y = y;
      obj.allEqual = size(obj.y, 1) == 0 || all(obj.y == obj.y(1, :));
      if obj.allEqual
        return
      end
      if ~isempty(obj.degree)
        obj.XP = generateFeatures(X, obj.degree, true);
        obj.minSize = max(obj.minSize, size(obj.XP, 2));
      end
      obj.current = obj.getData(true(n, 1));
    end
    
    function gain = get(obj, splitter)
    % evaluates splitter function
      if obj.allEqual
        gain = -inf;
        return;
      end
      idx = splitter(obj.X);
      if sum(idx) < obj.minSize || sum(~idx) < obj.minSize
        gain = -inf;
        return;
      end
      left = obj.getData(idx);
      right = obj.getData(~idx);
      n = size(obj.y, 1);
      nLeft = size(left.y, 1);
      nRight = size(right.y, 1);
      gain = obj.current.value ...
        - nLeft/n * left.value ...
        - nRight/n * right.value;
      % gain = obj.current.value - left.value - right.value;
    end
  end
  
  methods (Abstract, Access = protected)
    value = getValue(obj, data)
    % evaluates data using custom metric
  end
  
  methods (Access = protected)
    function data = getData(obj, idx)
    % returns portion of data given by idx and stores its value
      data = struct;
      data.idx = idx;
      data.y = obj.y(idx, :);
      if isempty(obj.modelFunc) && isempty(obj.degree)
        % constant model
        yPred = mean(data.y);
        sd2 = var(data.y);
        data.yPred = repmat(yPred, size(data.y));
        data.sd2 = repmat(sd2, size(data.y));
      elseif isempty(obj.modelFunc)
        % polynomial model
        X = obj.XP(idx, :);
        r = rank(X);
        % check rank deficiency
        if r < size(X, 2)
          d = size(obj.X, 2);
          if r >= d+1
            % use linear model
            X = obj.XP(idx, 1:d+1);
          else
            % use constant model
            X = obj.XP(idx, 1);
          end
        end
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:singularMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        M = (X' * X)^-1;
        b = M * X' * data.y;
        %b = X \ data.y;
        warning('on', 'MATLAB:rankDeficientMatrix');
        warning('on', 'MATLAB:singularMatrix');
        warning('on', 'MATLAB:nearlySingularMatrix');
        data.yPred = X * b;
        % var(b) = E(b^2) * (X'*X)^-1
        sd2 = immse(data.y, data.yPred);
        bCov = sd2 * M;
        data.sd2 = sum(X * bCov .* X, 2);
      else
        % custom model
        X = obj.X(idx, :);
        xMean = mean(X);
        model = obj.modelFunc(xMean);
        model = model.trainModel(X, data.y, xMean, 0);
        if obj.probabilistic
          [data.yPred, data.sd2, data.ci, data.p] = ...
            model.modelPredict(X, data.y);
        else
          [data.yPred, data.sd2] = ...
            model.modelPredict(X);
        end
      end
      data.value = obj.getValue(data);
    end
  end
end