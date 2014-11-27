function tests = gpModelTest
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../gpeda/src/vendor/gpml-matlab-v3.2/');
  run '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  
end

function test1DSinTest(testCase)
  X = (0:0.1:2)';
  y = sin(X);
  
  % Test this:
  m = GpModel([], 1);
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
  verifyLessThan(testCase, mae, 0.3);
  disp(['Gaussian process MAE = ' num2str(mae)]);
  verifyLessThan(testCase, dev, 0.5);
  disp(['Gaussian process mean std = ' num2str(mean(dev))]);
end

function testGetNTrainData(testCase)
  m = GpModel([], [0]);
  verifyEqual(testCase, m.getNTrainData, 3);
  m = GpModel([], [0 1 2 3 4]);
  verifyEqual(testCase, m.getNTrainData, 15);
end
