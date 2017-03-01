function tests = gpRfModelTest
  tests = functiontests(localfunctions);
end

function testGpRfModelConstructor(testCase)
  % test the GP+RF model constructor
  xMean = (0:0.1:2);
  
  % Test constructor works with empty settigs
  modelOpts = [];
  verifyNotEmpty(testCase, GpRfModel(modelOpts, xMean))
end

function testGpRfModelGetNTrainData(testCase)
  % test getNTrainData function
  dims = [1, 3, 5, 10, 20];
  nTrainData = 3*dims;
  
  for d = 1 : length(dims)
    xMean = ones(1, dims(d));
    
    % construct model
    modelOpts = [];
    m = GpRfModel(modelOpts, xMean);
    verifyEqual(testCase, m.getNTrainData, nTrainData(d))
  end
end

function testGpRfTrainModel(testCase)
% test gpRfModel train function for different dimensions
  dims = [2, 3, 5, 10, 20];
  generation = 1;
  
  % test different dimensions
  for d = 1:length(dims)
    rng(1);
    xTrain = randn(3*dims(d), dims(d));
    rng(2);
    yTrain = randn(3*dims(d), 1);
    xMean = mean(xTrain);
    
    % construct model
    modelOpts = [];
    m = GpRfModel(modelOpts, xMean);
    m = m.trainModel(xTrain, yTrain, xMean, generation);
    
    % trained model
    verifyNotEmpty(testCase, m)
  end

end

function testGpRfModelPredict(testCase)
% test gpRfModel train function for different dimensions
  dims = [2, 3, 5, 10, 20];
  generation = 1;
  fun = @(x) sin(sum(x, 2));
  
  % test different dimensions
  for d = 1:length(dims)
    rng(1);
    xTrain = randn(3*dims(d), dims(d));
    yTrain = fun(xTrain);
    xMean = mean(xTrain);
    rng(3);
    xTest = rand(dims(d)); 
    yTest = fun(xTest);
    
    % construct model
    modelOpts = [];
    m = GpRfModel(modelOpts, xMean);
    m = m.trainModel(xTrain, yTrain, xMean, generation);
    y = m.modelPredict(xTest);
    % trained model
    verifyGreaterThanOrEqual(testCase, y, -1)
    verifyLessThanOrEqual(testCase, y, 1)
  end

end