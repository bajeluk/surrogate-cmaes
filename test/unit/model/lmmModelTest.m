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
  m = LmmModel([], xMean);
  
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

function testKnn(testCase)
% test settings for knn

  dim = 2;
  nPoints = 20;
  f1 = @(x) sum(x.^2, 2);
  % add minimum to input
  X = randn(nPoints, dim);
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

  % test evaluation of knn option using char
  modelOptions.lmmKnn = 'obj.dim^3 + 1';
  % construct model
  m = LmmModel(modelOptions, xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % test k value
  verifyEqual(testCase, m.k, dim^3 + 1)

  % test evaluation of knn option using settings for all training points
  % Inf settings:
  modelOptions.lmmKnn = Inf;
  % construct model
  m = LmmModel(modelOptions, xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % test k value
  verifyEqual(testCase, m.k, nPoints)
  % NaN settings:
  modelOptions.lmmKnn = NaN;
  % construct model
  m = LmmModel(modelOptions, xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % test k value
  verifyEqual(testCase, m.k, nPoints)

  % test evaluation of knn option using incorrect settings
  modelOptions.lmmKnn = 0;
  % construct model
  m = LmmModel(modelOptions, xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % error is catched in Model, therefore the model should not be trained
  verifyFalse(testCase, m.isTrained())

end

function testGetTrainSet(testCase)
% test getTrainSet function
  dim = 2;
  X = 5 + randn(20, dim);
  y = fsphere(X);
  % nearest neighbors
  knnX = -1*ones(2, dim);
  knny = fsphere(knnX);
  X = [X; knnX];
  y = [y; knny];

  knn = 'obj.dim';
  xInput = randn(10, dim);
  sigma = 1;
  % following 3 lines are taken from CMA-ES code
  [B, diagD] = eig(cov(X));
  diagD = diag(sqrt(diagD)); % D contains standard deviations now
  BD = B*diag(diagD); % identical to BD = B.*repmat(diagD', dim, 1);
  xMean = mean(X);

  % set state variables
  stateVariables.xmean = xMean';
  stateVariables.sigma = sigma;
  stateVariables.lambda = 4 + floor(3*log(dim));
  stateVariables.BD = BD;
  stateVariables.diagD = diagD;
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

  % test
  modelOptions.lmmKnn = knn;
  % construct model
  m = LmmModel(modelOptions, xMean);
  % train model
  m = m.train(X, y, stateVariables, sampleOpts);
  % get training set
  [Xtr, ytr] = m.getTrainSet(xInput);
  % test same number of points
  verifyEqual(testCase, size(Xtr, 1), size(ytr, 1))

end

function y = fsphere(x)
  y = sum(x.^2, 2);
end