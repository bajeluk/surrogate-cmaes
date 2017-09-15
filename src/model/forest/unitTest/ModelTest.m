classdef (Abstract) ModelTest < Test
  
  properties (Constant, Access = private)
    resultsTemplate = struct(... % template for results
        'name', '', ...
        'params', struct, ...
        'trainRMSE', 0, ...
        'testRMSE', 0);
  end
  
  properties (Access = protected)
    results = [];
  end
  
  methods (TestMethodTeardown)
    function saveResults(testCase)   
      path = 'results';
      [~,~,~] = mkdir(path);
      filename = sprintf('%s/%s.mat', ...
        path, ...
        testCase.name{1});
      try
        load(filename);
        results = [results; testCase.results];
      catch
        results = testCase.results;
      end
      numel(results)
      save(filename, 'results');
    end
  end
  
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
      tic;
      model = model.trainModel(train.XN, train.y);
      time = toc;
      [train.yPred, train.sd2, train.ci] = ...
        model.modelPredict(train.XN);
      [test.yPred, test.sd2, test.ci] = ...
        model.modelPredict(test.XN);
      train.err = sqrt(immse(train.y, train.yPred));
      test.err = sqrt(immse(test.y, test.yPred));
      
      result = struct;
      result.name = sprintf('%s(%s)', ...
          testCase.name{2}, ...
          testCase.joinedParams);
      result.params = testCase.params;
      result.trainRMSE = train.err;
      result.testRMSE = test.err;
      if isempty(testCase.results)
        testCase.results = [result];
      else
        testCase.results(end+1) = result;
      end
      
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
    
    function [model, train, test, time] = testCoco(testCase, modelFunc, fNum, m)
      if nargin < 4
        m = 5;
      end
      % random points
      d = 2;
      n = 250 * d;
      X = testCase.generateInput(n, d, -m, m);
      y = benchmarks(X', fNum)';
      
      [model, train, test, time] = testCase.testModel(X, y, modelFunc);
    end
  end
end

