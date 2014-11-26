function tests = rfModelTest
  tests = functiontests(localfunctions);
end

function test1DSinTest(testCase)
  X = (0:0.1:2)';
  y = sin(X);
  
  % Test this:
  m = RfModel([], 1);
  m = m.train(X, y, 1, 1);
  [yPred, dev] = m.predict(X);
  
  mse = mean((y - yPred).^2);
  verifyLessThan(testCase, mse, 0.2);
  verifyLessThan(testCase, dev, 0.1);
end

function testGetNTrainData(testCase)
  m = RfModel([], [0]);
  verifyEqual(testCase, m.getNTrainData, 5);
  m = RfModel([], [0 1 2 3 4]);
  verifyEqual(testCase, m.getNTrainData, 25);
end
