classdef TreeModelTest < ModelTest
  
  properties (TestParameter)
    testMethod = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
  end
  
  methods (Test)
    function testConstantFunctionConstantPredictor(testCase)
      params = struct;
      testCase.reset(params);
      
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
    end
    
    function testLinearFunctionConstantPredictor(testCase)
      params = struct;
      testCase.reset(params);
      
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
    end
    
    function testLinearFunctionConstantPredictorPcaSplit(testCase)
      params = struct;
      testCase.reset(params);
      
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      
      modelOptions = struct;
      modelOptions.splits = {...
        AxisSplit(struct('transformationOptions', struct('pca', true)))...
        };
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testConstantFunctionLinearPredictor(testCase)
      params = struct;
      testCase.reset(params);
      
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      predictorOptions = struct('modelSpec', 'linear');
      predictorFunc = @() PolynomialModel(predictorOptions);
      splitGainOptions = struct('degree', 'linear');
      splitGain = SSESplitGain(splitGainOptions);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitGain', splitGain ...
        );
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testLinearFunctionLinearPredictor(testCase)
      params = struct;
      testCase.reset(params);
      
      % random points
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      % two linear functions
      split = X(:, 1) <= m/2;
      y = [5 + 2*X(:, 1) + 3*X(:, 2)] .* split ...
        + [5 + 3*X(:, 1) + 2*X(:, 2)] .* ~split;
      
      predictorOptions = struct('modelSpec', 'linear');
      predictorFunc = @() PolynomialModel(predictorOptions);
      splitGainOptions = struct('degree', 'linear');
      splitGain = SSESplitGain(splitGainOptions);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitGain', splitGain ...
        );
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testQuadraticFunctionQuadraticPredictor(testCase)
      params = struct;
      testCase.reset(params);
      
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
      splitGainOptions = struct('degree', 'quadratic');
      splitGain = SSESplitGain(splitGainOptions);
      modelOptions = struct(...
        'predictorFunc', predictorFunc, ...
        'splitGain', splitGain ...
        );
      modelFunc = @() TreeModel(modelOptions);
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
  end
end