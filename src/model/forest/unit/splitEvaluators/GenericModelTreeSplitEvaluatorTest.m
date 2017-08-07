classdef GenericModelTreeSplitEvaluatorTest < matlab.unittest.TestCase
  methods (Test)
    function testConstantGain(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      TreeSplitEvaluatorTest.draw(X, y, split);
      
      modelOptions = struct;
      modelFunc = @(xMean) ConstantModel(modelOptions, xMean);
      evaluator = GenericModelTreeSplitEvaluator(@gainMse, modelFunc);
      evaluator.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/4;
      gain = evaluator.eval(splitter);
      % small gain
      verifyLessThanOrEqual(testCase, gain, 0.1);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m*3/4;
      gain = evaluator.eval(splitter);
      % small gain
      verifyLessThanOrEqual(testCase, gain, 0.1);
      
      % try to split in half in x axis
      splitter = @(X) X(:, 1) <= m/2;
      gain = evaluator.eval(splitter);
      % bigger gain
      verifyGreaterThanOrEqual(testCase, gain, 0.24);
    end
    
    function testLinearNoGain(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      TreeSplitEvaluatorTest.draw(X, y, true(n, 1));
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      evaluator = GenericModelTreeSplitEvaluator(@gainMse, modelFunc);
      evaluator.reset(X, y);
      
      % try to split in half in x axis
      splitter = @(X) X(:, 1) <= m/2;
      gain = evaluator.eval(splitter);
      
      % there's no gain in splitting it
      verifyLessThanOrEqual(testCase, gain, 1e-10);
    end
    
    function testLinearGain(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      TreeSplitEvaluatorTest.draw(X, y, split);
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      evaluator = GenericModelTreeSplitEvaluator(@gainMse, modelFunc);
      evaluator.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/4;
      gain = evaluator.eval(splitter);
      % some gain is there
      verifyGreaterThanOrEqual(testCase, gain, 100);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m*3/4;
      gain = evaluator.eval(splitter);
      % some gain is there
      verifyGreaterThanOrEqual(testCase, gain, 90);
      
      % try to split in half in x axis
      splitter = @(X) X(:, 1) <= m/2;
      gain = evaluator.eval(splitter);
      % bigger gain
      verifyGreaterThanOrEqual(testCase, gain, 250);
    end
    
    function testLinearGainOnQuadraticFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two quadratic functions
      split = X(:, 1) <= m/2;
      y = [2*X(:, 1).^2 + 3*X(:, 2).^2] .* split ...
        + [3*X(:, 1).^2 + 2*X(:, 2).^2] .* ~split;
      TreeSplitEvaluatorTest.draw(X, y, split);
      
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      evaluator = GenericModelTreeSplitEvaluator(@gainMse, modelFunc);
      evaluator.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/4;
      gain = evaluator.eval(splitter);
      % huge gain
      verifyGreaterThanOrEqual(testCase, gain, 4e6);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m*3/4;
      gain = evaluator.eval(splitter);
      % huge gain
      verifyGreaterThanOrEqual(testCase, gain, 5e6);
      
      % try to split in half in x axis
      splitter = @(X) X(:, 1) <= m/2;
      gain = evaluator.eval(splitter);
      % enormous gain
      verifyGreaterThanOrEqual(testCase, gain, 7e6);
    end
    
    function testQuadraticGainOnQuadraticFunction(testCase)
      rng('default');
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two quadratic functions
      split = X(:, 1) <= m/2;
      y = [2*X(:, 1).^2 + 3*X(:, 2).^2] .* split ...
        + [3*X(:, 1).^2 + 2*X(:, 2).^2] .* ~split;
      TreeSplitEvaluatorTest.draw(X, y, split);
      
      modelOptions = struct('modelSpec', 'quadratic');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      evaluator = GenericModelTreeSplitEvaluator(@gainMse, modelFunc);
      evaluator.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/4;
      gain = evaluator.eval(splitter);
      % big gain
      verifyGreaterThanOrEqual(testCase, gain, 2e5);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m*3/4;
      gain = evaluator.eval(splitter);
      % big gain
      verifyGreaterThanOrEqual(testCase, gain, 1e5);
      
      % try to split in half in x axis
      splitter = @(X) X(:, 1) <= m/2;
      gain = evaluator.eval(splitter);
      % enormous gain
      verifyGreaterThanOrEqual(testCase, gain, 8e5);
    end
  end
end