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
  m = m.train(X, y, 1, 1);
  [yPred, dev] = m.predict(X);
  
  mse = mean((y - yPred).^2);
  verifyLessThan(testCase, mse, 0.2);
  verifyLessThan(testCase, dev, 0.1);
end

function testGetNTrainData(testCase)
  m = GpModel([], [0]);
  verifyEqual(testCase, m.getNTrainData, 3);
  m = GpModel([], [0 1 2 3 4]);
  verifyEqual(testCase, m.getNTrainData, 15);
end
