classdef TreeModelTest < ModelTest
  
  methods (Test)
    function testConstantFunctionConstantPredictor(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      modelOptions = struct;
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
    
    function testLinearFunctionConstantPredictor(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      
      modelOptions = struct;
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % small error
      verifyLessThan(testCase, train.err, 50);
      verifyLessThan(testCase, test.err, 160);
    end
    
    function testLinearFunctionConstantPredictorPcaSplit(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      
      modelOptions = struct;
      modelOptions.splitGenerators = {...
        AxisTreeSplitGenerator(struct('transformationOptions', struct('pca', true)))...
        };
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % better with pca
      verifyLessThan(testCase, train.err, 20);
      verifyLessThan(testCase, test.err, 50);
    end
    
    function testConstantFunctionLinearPredictor(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      predictorOptions = struct('modelSpec', 'linear');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitEvaluator', splitEvaluator ...
        );
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
    
    function testLinearFunctionLinearPredictor(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      
      predictorOptions = struct('modelSpec', 'linear');
      predictorFunc = @(xMean) PolynomialModel(predictorOptions, xMean);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitEvaluator', splitEvaluator ...
        );
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
    
    function testQuadraticFunctionQuadraticPredictor(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two quadratic functions
      split = X(:, 1) <= m/2;
      y = [2*X(:, 1).^2 + 3*X(:, 2).^2] .* split ...
        + [3*X(:, 1).^2 + 2*X(:, 2).^2] .* ~split;
      
      predictorOptions = struct('modelSpec', 'quadratic');
      predictorFunc = @() PolynomialModel(predictorOptions);
      splitEvaluator = PolynomialModelTreeSplitEvaluator(...
        @gainMse, predictorOptions.modelSpec);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitEvaluator', splitEvaluator ...
        );
      modelFunc = @(xMean) TreeModel(modelOptions, xMean);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
  end
end