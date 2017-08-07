classdef ConstantModelTest < matlab.unittest.TestCase
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
      
      modelOptions = struct;
      modelFunc = @(xMean) ConstantModel(modelOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      % some unavoidable error
      verifyLessThan(testCase, abs(train.err - 0.25), 1e-2);
      verifyLessThan(testCase, abs(test.err - 0.25), 1e-2);
    end
    
    function testConstantFunctionSet(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      modelOptions = struct('y', 0.75);
      modelFunc = @(xMean) ConstantModel(modelOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      % bigger error
      verifyLessThan(testCase, abs(train.err - 0.30), 1e-2);
      verifyLessThan(testCase, abs(test.err - 0.31), 1e-2);
    end
    
    function testLinearFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      
      modelOptions = struct;
      modelFunc = @(xMean) ConstantModel(modelOptions, xMean);
      [model, train, test, time] = ModelTest.testModel(X, y, modelFunc);
      
      % big error
      verifyGreaterThan(testCase, train.err, 1e4);
      verifyGreaterThan(testCase, test.err, 1e4);
    end
  end
end