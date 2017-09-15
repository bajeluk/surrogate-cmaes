classdef TreeModelTestFuzziness < ModelTest
  
  properties (TestParameter)
    fNum = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    m = {10}; %{5, 10, 100};
    modelSpec = {'linear', 'quadratic'};
    minLeafSize = {10};
    minGain = {0.1};
    splitGain1 = {'SSE'};
    splitGain2 = {'SSE'};
    fuzziness = {0, 0.01, 0.1, 0.25, 0.4};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = true;
    end
  end
  
  methods (Test)
    function test1(testCase, fNum, m, ...
        minLeafSize, minGain, splitGain1, fuzziness)
      params = struct;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain1;
      params.fuzziness = fuzziness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
      splitGain = splitGainFunc(splitGainOptions);
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = 1 / fuzziness;
      splits{1} = AxisSplit(splitOptions);
      
      modelOptions = struct;
      modelOptions.minGain = minGain;
      modelOptions.splitGain = splitGain;
      modelOptions.splits = splits;
      modelOptions.fuzziness = fuzziness;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
    
    function test2(testCase, fNum, m, ...
        modelSpec, minLeafSize, minGain, splitGain2, fuzziness)
      params = struct;
      params.modelSpec = modelSpec;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain2;
      params.fuzziness = fuzziness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      predictorOptions = struct;
      predictorOptions.modelSpec = modelSpec;
      predictorFunc = @() PolynomialModel(predictorOptions);
      
      splitGainOptions = struct;
      splitGainOptions.degree = modelSpec;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain2));
      splitGain = splitGainFunc(splitGainOptions);
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = 1 / fuzziness;
      splits{1} = AxisSplit(splitOptions);
      
      modelOptions = struct;
      modelOptions.predictorFunc = predictorFunc;
      modelOptions.minGain = minGain;
      modelOptions.splitGain = splitGain;
      modelOptions.splits = splits;
      modelOptions.fuzziness = fuzziness;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end