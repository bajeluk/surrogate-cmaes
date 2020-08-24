function tests = randomModelTest
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
% test random model on empty input

  % construct model with empty options
  modelOpts = [];
  m = RandomModel(modelOpts, 1);
  % train model using empty input should raise warning
  trainHandle = @() m.train([], []);
  verifyWarning(testCase, trainHandle, 'scmaes:model:emptytrainset')
end

function test5DSinTest(testCase)
  % test the GP model prediction on 5D sin(x) function
  rng(5)
  X = 3*rand(30,5) - 1;
  y = sin(sqrt(sum(X.^2,2)));
  
  % Train on the train data:
  m = RandomModel([], mean(X,1));
  m = m.train(X, y);

  % Predict the test data:
  Xtest = 3*rand(100,5) - 1;
  ytest = sin(sqrt(sum(Xtest.^2,2)));
  [yPred, dev] = m.predict(Xtest);
  
  % verify size
  verifySize(testCase, yPred, size(ytest))
  verifySize(testCase, dev, size(ytest))
  % verify that data are from the range of y
  verifyTrue(testCase, max(yPred)<=max(y) && min(yPred)>=min(y))
end

function testGetNTrainData(testCase)
  % test of the helper function getNTrainData
  m = RandomModel([], 0);
  verifyEqual(testCase, m.getNTrainData, 2);
  m = RandomModel([], [0 1 2 3 4]);
  verifyEqual(testCase, m.getNTrainData, 10);
end