function tests = gpModelTest
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../gpeda/src/vendor/gpml-matlab-v3.2/');
  run '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  
end

function test1DSinTest(testCase)
  % test the GP model prediction on 1D sin(x) function
  X = (0:0.1:2)';
  y = sin(X);
  
  % Test this:
  modelOpts = [];
  m = GpModel(modelOpts, 1);
  warning('off');
  m = m.train(X, y, 1, 1);
  warning('on');
  [yPred, dev] = m.predict(X);
  
  mse = mean((y - yPred).^2);
  verifyLessThan(testCase, mse, 0.2);
  disp(['Gaussian process MSE = ' num2str(mse)]);
  verifyLessThan(testCase, dev, 0.1);
  disp(['Gaussian process mean std = ' num2str(mean(dev))]);
end

function test5DSinTest(testCase)
  % test the GP model prediction on 5D sin(x) function
  X = 3*rand(30,5) - 1;
  y = sin(sqrt(sum(X.^2,2)));
  
  % Train on the train data:
  m = GpModel([], mean(X,1));
  warning('off');
  m = m.train(X, y, 1, 1);
  warning('on');

  % Train on the test data:
  Xtest = 3*rand(100,5) - 1;
  ytest = sin(sqrt(sum(Xtest.^2,2)));
  [yPred, dev] = m.predict(Xtest);
  
  mse = mean((ytest - yPred).^2);
  verifyLessThan(testCase, mse, 0.3);
  disp(['Gaussian process MSE = ' num2str(mse)]);
  mae = mean(abs(ytest - yPred));
  verifyLessThan(testCase, mae, 0.35);
  disp(['Gaussian process MAE = ' num2str(mae)]);
  verifyLessThan(testCase, dev, 0.5);
  disp(['Gaussian process mean std = ' num2str(mean(dev))]);
end

function testGetNTrainData(testCase)
  % test of the helper function getNTrainData
  m = GpModel([], [0]);
  verifyEqual(testCase, m.getNTrainData, 3);
  m = GpModel([], [0 1 2 3 4]);
  verifyEqual(testCase, m.getNTrainData, 15);
end

function test1DLinearShift(testCase)
  % test of the shift() function for Generation-based EC
  % without original evaluations for aligning the model
  X = (0:0.5:2)';
  y = linear(X);
  
  % Train the model
  m = GpModel([], 1);
  warning('off');
  m = m.train(X, y, 1, 1);
  warning('on');
  
  % Shift the new mean on the boundary of the dataset
  shift = [1];
  m = m.shift(mean(X,1)+shift);
  xShift = X + repmat(shift,size(X,1),1);
  yShift = linear(xShift);

  % Predict with the shifted model
  [yPred, dev] = m.predict(xShift);

  mse = mean((yShift - yPred).^2);
  verifyLessThan(testCase, mse, 0.2);
  disp(['Gaussian process MSE = ' num2str(mse)]);
  verifyLessThan(testCase, dev, 0.1);
  disp(['Gaussian process mean std = ' num2str(mean(dev))]);
end

function test1DLinearShiftReevaluate(testCase)
  % test of the shiftReevaluate() function for Generation-based EC
  % on the 1D linear function
  % WITH the original evaluations for aligning the model
  X = (0:0.5:2)';
  y = linear(X);
  
  % Train the model
  m = GpModel([], 1);
  warning('off');
  m = m.train(X, y, 1, 1);
  warning('on');
  
  % Shift the new mean -- this can be far away! :)
  shift = [2];
  [m, evals] = m.shiftReevaluate(mean(X,1)+shift, @linear);
  verifyEqual(testCase, evals, 1);
  xShift = X + repmat(shift,size(X,1),1);
  yShift = linear(xShift);

  % Predict with the shifted model
  [yPred, dev] = m.predict(xShift);

  mse = mean((yShift - yPred).^2);
  verifyLessThan(testCase, mse, 0.2);
  disp(['Gaussian process MSE = ' num2str(mse)]);
  verifyLessThan(testCase, dev, 0.1);
  disp(['Gaussian process mean std = ' num2str(mean(dev))]);
end

function testShiftReevaluateNaNs(testCase)
  % test if the shiftReevaluate() function handles NaNs from
  % the fitness correctly
  X = (0:0.5:2)';
  y = linear(X);
  
  % Train the model
  m = GpModel([], 1);
  warning('off');
  m = m.train(X, y, 1, 1);
  warning('on');
  
  % Try if it returns zero evaluations when fitness
  % returns all-the-times NaN
  shift = [2];
  [m, evals] = m.shiftReevaluate(mean(X,1)+shift, @nanFunction);
  verifyEqual(testCase, evals, 0);

  % Test if fitness which returns once NaN works
  global nans;
  nans = 0;
  [~, evals] = m.shiftReevaluate(mean(X,1)+shift, @getOnceNaN);
  verifyEqual(testCase, evals, 1);
end


% Simply fitness functions for this test

function y = linear(x)
  y = -0.5 * x + 5;
end

function y = nanFunction(x)
  y = nan(size(x,1),1);
end

function y = getOnceNaN(x)
  global nans;
  if (nans == 0)
    y = nan(size(x,1),1);
    nans = 1;
  else
    y = ones(size(x,1),1);
  end
end
