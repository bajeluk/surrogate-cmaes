classdef ModelTest
  properties (Constant)
    drawEnabled = true
  end
  
  methods (Static)
    function [model, train, test, time] = testModel(X, y, modelFunc)
      n = size(X, 1);
      idx = rand(n, 1) < 3/4;
      train = struct('X', X(idx, :), 'y', y(idx, :));
      test = struct('X', X(~idx, :), 'y', y(~idx, :));
      normalizeFunc = @(X1) X1; %...
        %bsxfun(@rdivide, bsxfun(@minus, X1, mean(train.X)), std(train.X));
      train.XN = normalizeFunc(train.X);
      test.XN = normalizeFunc(test.X);
      
      xMean = mean(train.XN);
      model = modelFunc(xMean);
      tic
      model.trainModel(train.XN, train.y, xMean, 0);
      time = toc;
      [train.yPred, train.sd2] = model.modelPredict(train.XN);
      [test.yPred, test.sd2] = model.modelPredict(test.XN);
      train.err = immse(train.y, train.yPred);
      test.err = immse(test.y, test.yPred);
      
      if ModelTest.drawEnabled
        stack = dbstack;
        figure('Name', stack(2).name);
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
        legend('y', 'y_{pred}^{train}', 'y_{pred}^{test}');
      end
    end
  end
end

