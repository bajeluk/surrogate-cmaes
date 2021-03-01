function tests = lmmModelTest
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
  % test the local meta-model with empty input
  modelOpts = [];
  X = [];
  y = [];
  
  m = LmmModel(modelOpts, 1);
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
  % set state variables
  stateVariables.xmean = xMean';
  stateVariables.sigma = 1;
  stateVariables.lambda = 4 + floor(3*log(dim));
  stateVariables.BD = eye(dim);
  stateVariables.diagD = ones(1, dim);
  stateVariables.countiter = 1;
  % set sampleOpts
  sampleOpts.noiseReevals = 0;
  sampleOpts.lbounds = -5;
  sampleOpts.ubounds = 5;
  sampleOpts.isBoundActive = true;
  sampleOpts.counteval = 1;
  sampleOpts.flgEvalParallel = false;
  sampleOpts.flgDiagonalOnly = false;
  sampleOpts.noiseHandling = 0;
  sampleOpts.xintobounds = @xintobounds;
  
  % construct model
  m = LmmModel([], xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % get model minimum
  mMin = m.minimumX();
  
  verifyEqual(testCase, mMin, xmin, 'AbsTol', 1e-8)
  
end