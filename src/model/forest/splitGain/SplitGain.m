classdef SplitGain
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
  
  properties %(Access = protected)
    X % input data
    y % output data
    evalFunc % function which evaluates data
    current % struct holding data for current node
    modelFunc % function which creates new model, model = modelFunc(xMean)
    probabilistic % whether to compute also probability and confidence interval
    XP % generated polynomial features for MethodPolynomialModel
    degree % degree of polynomial features for MethodPolynomialModel
  end
  
  methods
    function obj = SplitGain(evalFunc, modelSpec, probabilistic)
    % constructor
    % evaluator = SplitGain()
    %   evaluates functions based on the data
    %   (y, variance)
    % evaluator = SplitGain(degree)
    %   evaluates functions based on the output of polynomial model with
    %   given degree
    %   (y, yPred, variancePred)
    % evaluator = SplitGain(modelFunc)
    %   evaluates functions based on the output of the specified model
    %   (y, yPred, variancePred)
    % evaluator = SplitGain(modelFunc, true)
    %   evaluates functions based on the output of the specified model
    %   (y, yPred, variancePred, confidenceIntervalPred, probabilityPred)
      if nargin >= 1
        obj.evalFunc = evalFunc;
      end
      if nargin >= 2
        if isa(modelSpec, 'function_handle')
          obj.modelFunc = modelSpec;
        else
          obj.degree = modelSpec;
        end
      end
      if nargin >= 3
        obj.probabilistic = probabilistic;
      end
    end
    
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      [n, ~] = size(X);
      obj.X = X;
      obj.y = y;
      if ~isempty(obj.degree)
        obj.XP = generateFeatures(X, obj.degree, true);
      end
      obj.current = obj.getData(true(n, 1));
    end
    
    function gain = eval(obj, splitter)
    % evaluates splitter function
      idx = splitter(obj.X);
      left = obj.getData(idx);
      right = obj.getData(~idx);
      gain = obj.current.value - left.value - right.value;
    end
  end
  
  methods (Access = protected)
    function value = getValue(obj, data)
    % evaluates data using custom metric
      if isempty(obj.evalFunc)
        value = 0;
      else
        value = obj.evalFunc(data);
      end
    end
    
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
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:singularMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        data.yPred = X * (X \ data.y);
        warning('on', 'MATLAB:rankDeficientMatrix');
        warning('on', 'MATLAB:singularMatrix');
        warning('on', 'MATLAB:nearlySingularMatrix');
        sd2 = mean((data.y - data.yPred).^2);
        data.sd2 = repmat(sd2, size(data.y));
      else
        % custom model
        X = obj.X(idx, :);
        xMean = mean(X);
        model = obj.modelFunc(xMean);
        model.trainModel(X, data.y, xMean, 0);
        if obj.probabilistic
          [data.yPred, data.sd2, data.ci, data.p] = ...
            model.predict(X, data.y);
        else
          [data.yPred, data.sd2] = ...
            model.predict(X);
        end
      end
      data.value = obj.getValue(data);
    end
  end
end