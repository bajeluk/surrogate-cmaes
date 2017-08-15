classdef (Abstract) SplitGain
% SplitGain evaluates split functions used in decision trees using
% custom metric. The evaluation can be based on the data alone or on the
% output of a model. If the model is polynomial, optimized evaluation is
% used.
    
  enumeration
    MethodData, MethodModel, MethodProbabilityModel, MethodPolynomialModel
  end
  
  properties %(Access = protected)
    X % input data
    y % output data
    current % struct holding data for current node
    method = MethodData % method of evaluation
    modelFunc % function which creates new model, model = modelFunc(xMean)
    XP % generated polynomial features for MethodPolynomialModel
    degree % degree of polynomial features for MethodPolynomialModel
  end
  
  methods
    function obj = reset(obj, X, y)
    % resets the evaluator with new data (X, y)
      [n, ~] = size(X);
      obj.X = X;
      obj.y = y;
      switch obj.method
        case MethodPolynomialModel
          obj.XP = generateFeatures(X, obj.degree, true);
      end
      obj.current = obj.getData(true(n, 1));
    end
    
    function gain = eval(obj, splitter)
    % evaluates splitter function
      idx = splitter(obj.X);
      left = getData(idx);
      right = getData(~idx);
      gain = obj.current.value - (left.value + right.value);
    end
  end
  
  methods (Abstract, Access = protected)
    value = getValue(data)
    % evaluates data using custom metric
  end
  
  methods (Access = protected)
    function obj = SplitGain(modelSpec, probabilistic)
    % constructor
    % evaluator = SplitGain()
    %   evaluates functions based on the data
    %   (y, variance)
    % evaluator = SplitGain(modelFunc)
    %   evaluates functions based on the output of the specified model
    %   (y, yPred, variancePred)
    % evaluator = SplitGain(modelFunc, true)
    %   evaluates functions based on the output of the specified model
    %   (y, yPred, variancePred, confidenceIntervalPred, probabilityPred) 
      if nargin == 0
        obj.method = SplitGain.MethodData;
      elseif nargin >= 1
        if isa(modelSpec, 'function_handle')
          obj.modelFunc = modelSpec;
          if nargin >= 2 && probabilistic
            obj.method = SplitGain.MethodProbabilityModel;
          else
            obj.method = SplitGain.MethodModel;
          end
        else
          obj.method = SplitGain.MethodPolynomialModel;
          obj.degree = modelSpec;
        end
      end
    end
    
    function data = getData(obj, idx)
    % returns portion of data given by idx and stores its value
      data = struct;
      data.idx = idx;
      data.y = obj.y(idx, :);
      switch obj.method
        case SplitGain.MethodData
          sd2 = var(data.y);
          data.sd2 = repmat(sd2, size(data.y));
        case SplitGain.MethodPolynomialModel
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
        case SplitGain.MethodModel
          X = obj.X(idx, :);
          xMean = mean(X);
          model = obj.modelFunc(xMean);
          model.trainModel(X, data.y, xMean, 0);
          [data.yPred, data.sd2] = ...
            model.predict(X);
        case SplitGain.MethodProbabilityModel
          X = obj.X(idx, :);
          xMean = mean(X);
          model = obj.modelFunc(xMean);
          model.trainModel(X, data.y, xMean, 0);
          [data.yPred, data.sd2, data.ci, data.p] = ...
            model.predict(X, data.y);
      end
      data.value = obj.getValue(data);
    end
  end
end