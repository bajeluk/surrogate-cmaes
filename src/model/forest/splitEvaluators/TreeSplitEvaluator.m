classdef (Abstract) TreeSplitEvaluator < handle
  enumeration
    MethodData, MethodModel, MethodProbabilityModel, MethodPolynomialModel
  end
  
  properties (Access = protected)
    X % input data
    y % output data
    current % struct holding data in current node
    method
    modelFunc
    XP
    degree
    best
  end
  
  methods
    function reset(obj, X, y)
      % resets the evaluator with new data
      obj.X = X;
      obj.y = y;
      switch obj.method
        case MethodPolynomialModel
          obj.XP = generateFeatures(X, obj.degree, true);
      end
      obj.best = struct('gain', -inf);
      obj.resetData();
      obj.current = obj.getData(true(size(X, 1), 1));
    end
    
    function gain = eval(obj, f)
      % evaluates the splitting function
      idx = f(obj.X);
      left = getData(idx);
      right = getData(~idx);
      gain = obj.current.value - (left.value + right.value);
      if gain > obj.best.gain
        obj.best.gain = gain;
        obj.best.splitter = f;
      end
    end
    
    function best = get.best(obj)
      best = obj.best;
    end
  end
  
  methods (Abstract, Access = protected)
    value = getValue(data)
  end
  
  methods (Access = protected)
    function obj = TreeSplitEvaluator(varargin)
      if nargin == 1
        obj.method = TreeSplitEvaluator.MethodData;
      elseif nargin >= 2
        if isa(model, 'function_handle')
          if nargin >= 3 && varargin{3}
            obj.method = TreeSplitEvaluator.MethodProbabilityModel;
          else
            obj.method = TreeSplitEvaluator.MethodModel;
          end
          obj.modelFunc = varargin{1};
        else
          obj.method = TreeSplitEvaluator.MethodPolynomialModel;
          obj.degree = varargin{1};
        end
      end
    end
    
    function resetData(obj)
    end
    
    function data = getData(obj, idx)
      data = struct;
      data.idx = 1:numel(idx);
      data.idx = dta.idx(idx);
      data.y = obj.y(idx, :);
      switch obj.method
        case MethodData
          sd2 = var(data.y);
          data.sd2 = repmat(sd2, size(data.y));
        case MethodPolynomialModel
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
        case MethodModel
          X = obj.X(idx, :);
          xMean = mean(X);
          model = obj.modelFunc(xMean);
          model.trainModel(X, data.y, xMean, 0);
          [data.yPred, data.sd2] = ...
            model.predict(X);
        case MethodProbabilityModel
          X = obj.X(idx, :);
          xMean = mean(X);
          model = obj.modelFunc(xMean);
          model.trainModel(X, data.y, xMean, 0);
          [data.yPred, data.sd2, data.ci, data.p] = ...
            model.predict(X, data.y);
      end
      data.value = getValue(obj, data);
    end
  end
end