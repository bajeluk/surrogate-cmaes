classdef PolynomialModelTreeSplitEvaluator < TreeSplitEvaluator
  properties (Access = protected)
    degree % degree of polynomial
    XP % input data with polynomial features
  end
  
  methods
    function obj = PolynomialModelTreeSplitEvaluator(gainFunc, degree)
      if nargin > 0
        obj.gainFunc = gainFunc;
        obj.degree = degree;
      end
    end
    
    function reset(obj, X, y)
      % resets the evaluator with new data
      obj.X = X;
      obj.y = y;
      obj.current = struct;
      obj.current.y = y;
      obj.XP = generateFeatures(X, obj.degree, true);
      if isempty(obj.XP)
        obj.current.yPred = ...
          PolynomialModelTreeSplitEvaluator.modelPredict(obj.X, obj.y, obj.degree);
      else
        obj.current.yPred = ...
          PolynomialModelTreeSplitEvaluator.matrixPredict(obj.XP, obj.y);
      end
    end
    
    function r = eval(obj, f)
      idx = f(obj.X);
      
      left = struct;
      left.y = obj.y(idx);
      if isempty(obj.XP)
        X = obj.X(idx, :);
        left.yPred = ...
          PolynomialModelTreeSplitEvaluator.modelPredict(X, left.y, obj.degree);
      else
        X = obj.XP(idx, :);
        left.yPred = ...
          PolynomialModelTreeSplitEvaluator.matrixPredict(X, left.y);
      end
      
      right = struct;
      right.y = obj.y(~idx);
      if isempty(obj.XP)
        X = obj.X(~idx, :);
        right.yPred = ...
          PolynomialModelTreeSplitEvaluator.modelPredict(X, right.y, obj.degree);
      else
        X = obj.XP(~idx, :);
        right.yPred = ...
          PolynomialModelTreeSplitEvaluator.matrixPredict(X, right.y);
      end
      
      r = obj.gainFunc(obj.current, left, right);
    end
  end
  
  methods (Static, Access = private)
    function [yPred] = matrixPredict(X, y)
      warning('off', 'MATLAB:rankDeficientMatrix');
      warning('off', 'MATLAB:singularMatrix');
      warning('off', 'MATLAB:nearlySingularMatrix');
      yPred = X * (X \ y);
      warning('on', 'MATLAB:rankDeficientMatrix');
      warning('on', 'MATLAB:singularMatrix');
      warning('on', 'MATLAB:nearlySingularMatrix');
    end
    
    function [yPred] = modelPredict(X, y, degree)
      warning('off', 'stats:LinearModel:RankDefDesignMat');
      model = fitlm(X, y, degree);
      yPred = model.predict(X);
      warning('on', 'stats:LinearModel:RankDefDesignMat');
    end
  end
end