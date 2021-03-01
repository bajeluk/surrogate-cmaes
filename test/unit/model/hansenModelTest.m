function tests = hansenModelTest
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
  % test the linear-quadratic model with empty input
  modelOpts = [];
  X = [];
  y = [];
  
  m = HansenModel(modelOpts, 1);
  verifyNotEmpty(testCase, m);
  % empty training
  trainHandle = @() m.train(X, y);
  verifyWarning(testCase, trainHandle, 'scmaes:model:emptytrainset')
end

function testEmptySettings(testCase)
  % test default settings
  
  dim = 2;
  nPoints = 20;
  f1 = @(x) sum(x.^2, 2);
  X = randn(nPoints, dim);
  y = f1(X);
  xMean = mean(X);
  
  % construct model
  m = HansenModel([], xMean);
  
  verifyNotEmpty(testCase, m)
  
  % train model
  m = m.train(X, y);
  
  verifyTrue(testCase, m.isTrained())
end

function testMinimumX(testCase)
% test getting model minimum

  dim = 2;
  nPoints = 20;
  f1 = @(x) sum(x.^2, 2);
  xmin = zeros(1, dim);
  % add minimum to input
  X = [randn(nPoints-1, dim); xmin];
  y = f1(X);
  xMean = mean(X);
  
  % construct model
  m = HansenModel([], xMean);
  % train model
  m = m.train(X, y);
  % get model minimum
  mMin = m.minimumX();
  
  verifyEqual(testCase, mMin, xmin, 'AbsTol', 1e-14)
  
end