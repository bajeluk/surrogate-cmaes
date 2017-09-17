classdef TreeModelTestFuzziness < ModelTest
  
  properties (TestParameter)
    fNum = {1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    m = {10};
    modelSpec = {'linear', 'quadratic'};
    minLeafSize = {10};
    minGain = {0.1};
    splitGain0 = {'MSE'};
    splitGain1 = {'MSE'};
    %fuzziness = {0};
    %lambda = {0};
    fuzziness = {0.05, 1/6, 1/2};
    steepness = {1, 2, 5};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = true;
    end
  end
  
  methods (Test)
    function test0(testCase, fNum, m, ...
        minLeafSize, minGain, splitGain0, fuzziness, steepness)
      params = struct;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain0;
      params.fuzziness = fuzziness;
      params.steepness = steepness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain0));
      splitGain = splitGainFunc(splitGainOptions);
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = steepness;
      splits{1} = AxisSplit(splitOptions);
      
      modelOptions = struct;
      modelOptions.minGain = minGain;
      modelOptions.splitGain = splitGain;
      modelOptions.splits = splits;
      modelOptions.fuzziness = fuzziness;
      modelFunc = @() TreeModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
    
    function test1(testCase, fNum, m, ...
        modelSpec, minLeafSize, minGain, splitGain1, fuzziness, steepness)
      params = struct;
      params.modelSpec = modelSpec;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain1;
      params.fuzziness = fuzziness;
      params.steepness = steepness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      predictorOptions = struct;
      predictorOptions.modelSpec = modelSpec;
      predictorFunc = @() PolynomialModel(predictorOptions);
      
      splitGainOptions = struct;
      splitGainOptions.degree = modelSpec;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
      splitGain = splitGainFunc(splitGainOptions);
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = steepness;
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