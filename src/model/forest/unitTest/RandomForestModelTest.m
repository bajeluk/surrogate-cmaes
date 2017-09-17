classdef RandomForestModelTest < ModelTest

  properties (TestParameter)
    fNum = {2}; %{1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    m = {10};
    modelSpec = {'linear', 'quadratic'};
    minLeafSize = {10};
    minGain = {0.1};
    splitGain0 = {'MSE'};
    splitGain1 = {'MSE'};
    fuzziness = {0};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = true;
    end
  end
  
  methods (Test)
    function test0(testCase, fNum, m, ...
        minLeafSize, minGain, splitGain0, fuzziness)
      params = struct;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain0;
      params.fuzziness = fuzziness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain0));
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = 2;
      splits{1} = AxisSplit(splitOptions);
      
      treeModelOptions = struct;
      treeModelOptions.minGain = minGain;
      treeModelOptions.splitGain = splitGainFunc(splitGainOptions);
      treeModelOptions.splits = splits;
      treeModelOptions.fuzziness = fuzziness;
      treeModelFunc = @() TreeModel(treeModelOptions);
      
      modelOptions = struct;
      modelOptions.treeFunc = treeModelFunc;
      modelFunc = @() RandomForestModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
    
    function test(testCase, fNum, m, ...
        modelSpec, minLeafSize, minGain, splitGain1, fuzziness)
      params = struct;
      params.modelSpec = modelSpec;
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.splitGain = splitGain1;
      params.fuzziness = fuzziness;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      predictorOptions = struct;
      predictorOptions.modelSpec = modelSpec;
      predictorFunc = @() PolynomialModel(predictorOptions);
      
      splitGainOptions = struct;
      splitGainOptions.degree = modelSpec;
      splitGainOptions.minSize = minLeafSize;
      splitGainFunc = str2func(sprintf('%sSplitGain', splitGain1));
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = 2;
      splits{1} = AxisSplit(splitOptions);
      
      treeModelOptions = struct;
      treeModelOptions.predictorFunc = predictorFunc;
      treeModelOptions.minGain = minGain;
      treeModelOptions.splitGain = splitGainFunc(splitGainOptions);
      treeModelOptions.splits = splits;
      treeModelOptions.fuzziness = fuzziness;
      treeModelFunc = @() TreeModel(treeModelOptions);
      
      modelOptions = struct;
      modelOptions.treeFunc = treeModelFunc;
      modelFunc = @() RandomForestModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end