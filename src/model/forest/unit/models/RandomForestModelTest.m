classdef RandomForestModelTest < matlab.unittest.TestCase
  methods (Test)    
    function testSchwefelFunctionConstantPredictor(testCase)
      rng('default');
      % Schwefel function
      n = 1000;
      d = 2;
      f = @(A,B) 418.9829 * d - A .* sin(sqrt(abs(A))) - B .* sin(sqrt(abs(B)));
      X = rand(n, d) * 1000 - 500;
      y = f(X(:, 1), X(:, 2));
      
      treeOptions = struct;
      treeFunc = @(xMean) TreeModel(treeOptions, xMean);
      forestOptions = struct;
      forestOptions.treeFunc = treeFunc;
      modelFunc = @(xMean) RandomForestModel(forestOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      if ModelTest.drawEnabled
        [A, B] = meshgrid(-500:10:500, -500:10:500);
        C = f(A, B);
        hold on;
        surf(A, B, C);
        hold off;
      end
      
      % big error
      verifyLessThan(testCase, train.err, 1e4);
      verifyLessThan(testCase, test.err, 2e4);
    end
    
    function testSchwefelFunctionLinearPredictor(testCase)
      rng('default');
      % Schwefel function
      n = 1000;
      d = 2;
      f = @(A,B) 418.9829 * d - A .* sin(sqrt(abs(A))) - B .* sin(sqrt(abs(B)));
      X = rand(n, d) * 1000 - 500;
      y = f(X(:, 1), X(:, 2));
      
      predictorOptions = struct('modelSpec', 'linear');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      treeOptions = struct;
      treeOptions.minLeafSize = 5;
      treeOptions.predictorFunc = predictorFunc;
      treeOptions.splitEvaluator = splitEvaluator;
      treeFunc = @(xMean) TreeModel(treeOptions, xMean);
      forestOptions = struct;
      forestOptions.treeFunc = treeFunc;
      modelFunc = @(xMean) RandomForestModel(forestOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      if ModelTest.drawEnabled
        [A, B] = meshgrid(-500:10:500, -500:10:500);
        C = f(A, B);
        hold on;
        surf(A, B, C);
        hold off;
      end
      
      % big error
      verifyLessThan(testCase, train.err, 3e3);
      verifyLessThan(testCase, test.err, 4e4);
    end
    
    function testSchwefelFunctionQuadraticPredictor(testCase)
      rng('default');
      % Schwefel function
      n = 1000;
      d = 2;
      f = @(A,B) 418.9829 * d - A .* sin(sqrt(abs(A))) - B .* sin(sqrt(abs(B)));
      X = rand(n, d) * 1000 - 500;
      y = f(X(:, 1), X(:, 2));
      
      predictorOptions = struct('modelSpec', 'quadratic');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      treeOptions = struct;
      treeOptions.minLeafSize = 5;
      treeOptions.predictorFunc = predictorFunc;
      treeOptions.splitEvaluator = splitEvaluator;
      treeFunc = @(xMean) TreeModel(treeOptions, xMean);
      forestOptions = struct;
      forestOptions.treeFunc = treeFunc;
      modelFunc = @(xMean) RandomForestModel(forestOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      if ModelTest.drawEnabled
        [A, B] = meshgrid(-500:10:500, -500:10:500);
        C = f(A, B);
        hold on;
        surf(A, B, C);
        hold off;
      end
      
      % small train error, big test error
      verifyLessThan(testCase, train.err, 30);
      verifyGreaterThan(testCase, test.err, 3e4);
    end
    
    function testSchwefelFunctionQuadraticPredictorRbfSplitter(testCase)
      rng('default');
      % Schwefel function
      n = 1000;
      d = 2;
      f = @(A,B) 418.9829 * d - A .* sin(sqrt(abs(A))) - B .* sin(sqrt(abs(B)));
      X = rand(n, d) * 1000 - 500;
      y = f(X(:, 1), X(:, 2));
      
      predictorOptions = struct('modelSpec', 'quadratic');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      splitGenerator = RandomRbfTreeSplitGenerator('euclidean', 2, 10);
      treeOptions = struct;
      treeOptions.minLeafSize = 5;
      treeOptions.predictorFunc = predictorFunc;
      treeOptions.splitEvaluator = splitEvaluator;
      treeOptions.splitGenerators = {splitGenerator};
      treeFunc = @(xMean) TreeModel(treeOptions, xMean);
      forestOptions = struct;
      forestOptions.nTrees = 10;
      forestOptions.treeFunc = treeFunc;
      modelFunc = @(xMean) RandomForestModel(forestOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      if ModelTest.drawEnabled
        [A, B] = meshgrid(-500:10:500, -500:10:500);
        C = f(A, B);
        hold on;
        surf(A, B, C);
        hold off;
      end
      
      % small train error, big test error
      verifyLessThan(testCase, train.err, 30);
      verifyGreaterThan(testCase, test.err, 1e6);
    end
    
    function testSchwefelFunctionQuadraticPredictorSample(testCase)
      rng('default');
      % Schwefel function
      n = 1000;
      d = 2;
      f = @(A,B) 418.9829 * d - A .* sin(sqrt(abs(A))) - B .* sin(sqrt(abs(B)));
      X = rand(n, d) * 1000 - 500;
      y = f(X(:, 1), X(:, 2));
      
      predictorOptions = struct('modelSpec', 'quadratic');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      splitGenerator1 = RandomRbfTreeSplitGenerator('euclidean', 2, 10);
      splitGenerator2 = RandomPolynomialTreeSplitGenerator('quadratic', 2, 10);
      treeOptions = struct;
      treeOptions.minGain = 10;
      treeOptions.minLeafSize = 10;
      treeOptions.predictorFunc = predictorFunc;
      treeOptions.splitEvaluator = splitEvaluator;
      treeOptions.splitGenerators = {splitGenerator1, splitGenerator2};
      treeFunc = @(xMean) TreeModel(treeOptions, xMean);
      forestOptions = struct;
      forestOptions.nTrees = 10;
      forestOptions.treeFunc = treeFunc;
      modelFunc = @(xMean) RandomForestModel(forestOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      if ModelTest.drawEnabled
        [A, B] = meshgrid(-500:10:500, -500:10:500);
        C = f(A, B);
        hold on;
        surf(A, B, C);
        hold off;
      end
      
      % small error
      verifyLessThan(testCase, train.err, 7e2);
      verifyLessThan(testCase, test.err, 6e4);
    end
  end
end