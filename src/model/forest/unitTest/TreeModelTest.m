classdef TreeModelTest < ModelTest
  
  properties (TestParameter)
    %fNum = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    fNum = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, ...
      13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    m = {10}; %{5, 10, 100};
    modelSpec = {'linear', 'quadratic'};
    minLeafSize = {10, 50};
    minGain = {0.1, 1, 10};
    splitGain1 = {'DEMSD', 'DENN', 'DE', 'SSE', 'Var'};
    splitGain2 = {'DEMSD', 'DE', 'SSE', 'Var'};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = false;
    end
  end
  
  methods (Test)
    function test1(testCase, fNum, m, ...
        minLeafSize, minGain, splitGain1)
      params = struct;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain1;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
      splitGain = splitGainFunc(splitGainOptions);
      
      modelOptions = struct;
      modelOptions.minGain = minGain;
      modelOptions.splitGain = splitGain;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
    
    function test2(testCase, fNum, m, ...
        modelSpec, minLeafSize, minGain, splitGain2)
      params = struct;
      params.modelSpec = modelSpec;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain2;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      predictorOptions = struct;
      predictorOptions.modelSpec = modelSpec;
      predictorFunc = @() PolynomialModel(predictorOptions);
      
      splitGainOptions = struct;
      splitGainOptions.degree = modelSpec;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain2));
      splitGain = splitGainFunc(splitGainOptions);
      
      modelOptions = struct;
      modelOptions.predictorFunc = predictorFunc;
      modelOptions.minGain = minGain;
      modelOptions.splitGain = splitGain;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end