classdef (Abstract) WeakModelTest < ModelTest
  
  methods (Access = protected)
    function [model, train, test, time] = testModel(testCase, X, y, modelFunc)
      n = size(X, 1);
      idx = rand(n, 1) < 3/4;
      train = struct('X', X(idx, :), 'y', y(idx, :));
      test = struct('X', X(~idx, :), 'y', y(~idx, :));
      normalizeFunc = @(X1) X1; %...
        %bsxfun(@rdivide, bsxfun(@minus, X1, mean(train.X)), std(train.X));
      train.XN = normalizeFunc(train.X);
      test.XN = normalizeFunc(test.X);
      
      model = modelFunc();
      tic
      model = model.trainModel(train.XN, train.y);
      time = toc;
      [train.yPred, train.sd2, train.ci] = ...
        model.modelPredict(train.XN);
      [test.yPred, test.sd2, test.ci] = ...
        model.modelPredict(test.XN);
      train.err = sqrt(immse(train.y, train.yPred));
      test.err = sqrt(immse(test.y, test.yPred));
      
      if testCase.drawEnabled
        if size(X, 2) == 1
          plot(X, y);
          hold on;
          scatter(train.X, train.yPred);
          scatter(test.X, test.yPred);
          hold off;
        elseif size(X, 2) == 2
          scatter3(X(:, 1), X(:, 2), y);
          hold on;
          scatter3(train.X(:, 1), train.X(:, 2), train.yPred);
          scatter3(test.X(:, 1), test.X(:, 2), test.yPred);
          hold off;
        end
        if isempty(test.yPred)
          legend('y', 'y_{pred}^{train}');
        else
          legend('y', 'y_{pred}^{train}', 'y_{pred}^{test}');
        end
        description = sprintf('train RMSE: %.3f\n test RMSE: %.3f', ...
          train.err, test.err);
        title(description);
      end
    end
    
    function testConstantFunction(testCase, modelFunc)
      % random points
      d = 2;
      n = 250 * d;
      m = 100;
      X = rand(n, d) * m;
      % two constant functions
      split = X(:, 1) <= m/2;
      y = 1 * split;
      
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testLinearFunction(testCase, modelFunc)
      % random points
      d = 2;
      n = 250 * d;
      m = 100;
      X = rand(n, d) * m;;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);
      
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testQuadraticFunction(testCase, modelFunc)
      % random points
      d = 2;
      n = 250 * d;
      m = 100;
      X = rand(n, d) * m;
      % quadratic function
      y = 2*X(:, 1).^2 + 3*X(:, 2).^2;
      
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testDependentFeatures(testCase, modelFunc)
      % random points
      d = 2;
      n = 250 * d;
      m = 100;
      X = [ones(n, d-1), rand(n, 1) * m];
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);

      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
    
    function testFewPoints(testCase, modelFunc)
      % random points
      d = 2;
      n = 1 * d;
      m = 100;
      X = rand(n, d) * m;
      % linear function
      y = 5 + 2*X(:, 1) + 3*X(:, 2);

      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
  end
end

