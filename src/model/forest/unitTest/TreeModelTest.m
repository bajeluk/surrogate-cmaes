classdef TreeModelTest < ModelTest
  
  properties (TestParameter)
    fNum = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    %fNum = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, ...
    %  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    m = {5, 10, 100};
    modelSpec = {'constant', 'linear', 'quadratic'};
  end
  
  methods (Test)
    function test1(testCase, fNum, m)
      params = struct;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGain = SSESplitGain(splitGainOptions);
      
      modelOptions = struct;
      modelOptions.splitGain = splitGain;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
    
    function test(testCase, fNum, m, ...
        modelSpec)
      params = struct;
      params.modelSpec = modelSpec;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      predictorOptions = struct;
      predictorOptions.modelSpec = modelSpec;
      predictorFunc = @() PolynomialModel(predictorOptions);
      
      splitGainOptions = struct;
      splitGainOptions.degree = modelSpec;
      splitGain = SSESplitGain(splitGainOptions);
      
      modelOptions = struct;
      modelOptions.predictorFunc = predictorFunc;
      modelOptions.splitGain = splitGain;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end