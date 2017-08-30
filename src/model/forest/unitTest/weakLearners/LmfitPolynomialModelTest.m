classdef LmfitPolynomialModelTest < WeakModelTest
  methods (Test)
    function testConstantFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @() LmfitPolynomialModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % small error
      verifyLessThan(testCase, abs(train.err - 0.06), 1e-2);
      verifyLessThan(testCase, abs(test.err - 0.06), 1e-2);
    end
    
    function testLinearFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @() LmfitPolynomialModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
    
    function testQuadraticFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % quadratic function
      y = 2*X(:, 1).^2 + 3*X(:, 2).^2;
      
      modelOptions = struct('modelSpec', 'quadratic');
      modelFunc = @() LmfitPolynomialModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % no error
      verifyLessThan(testCase, train.err, 1e-4);
      verifyLessThan(testCase, test.err, 1e-4);
    end
    
    function testDependentFeatures(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = [ones(n, 1), rand(n, 1) * m];
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @() LmfitPolynomialModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % small error
      verifyLessThan(testCase, abs(train.err - 0), 1e-2);
      verifyLessThan(testCase, abs(test.err - 0), 1e-2);
    end
    
    function testFewPoints(testCase)
      rng('default');
      % random points
      n = 2;
      m = 100;
      X = rand(n, 2) * m;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @() LmfitPolynomialModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
      
      % small error
      verifyLessThan(testCase, abs(train.err - 0), 1e-2);
    end
  end
end