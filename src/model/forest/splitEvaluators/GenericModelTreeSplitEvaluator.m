classdef GenericModelTreeSplitEvaluator < TreeSplitEvaluator
  properties (Access = protected)
    modelFunc % function which creates a model (implements ModelImpl interface)
  end
  
  methods
    function obj = GenericModelTreeSplitEvaluator(gainFunc, modelFunc)
      if nargin > 0
        obj.gainFunc = gainFunc;
        obj.modelFunc = modelFunc;
      end
    end
    
    function reset(obj, X, y)
      % resets the evaluator with new data
      obj.X = X;
      obj.y = y;
      obj.current = struct;
      obj.current.y = y;
      xMean = mean(X);
      model = obj.modelFunc(xMean);
      model.trainModel(X, y, xMean, 0);
      obj.current.yPred = model.modelPredict(X);
    end
    
    function r = eval(obj, f)
      idx = f(obj.X);
      
      left = struct;
      left.y = obj.y(idx);
      X = obj.X(idx, :);
      xMean = mean(X);
      model = obj.modelFunc(xMean);
      model.trainModel(X, left.y, xMean, 0);
      left.yPred = model.modelPredict(X);
      
      right = struct;
      right.y = obj.y(~idx);
      X = obj.X(~idx, :);
      xMean = mean(X);
      model = obj.modelFunc(xMean);
      model.trainModel(X, right.y, xMean, 0);
      right.yPred = model.modelPredict(X);
      
      r = obj.gainFunc(obj.current, left, right);
    end
  end
end