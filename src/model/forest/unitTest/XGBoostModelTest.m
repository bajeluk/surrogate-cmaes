classdef XGBoostModelTest < ModelTest

  properties (TestParameter)
    fNum = {2}; %{1, 2, 6, 8, 13, 14, 15, 17, 20, 21};
    m = {10};
    minLeafSize = {10};
    minGain = {0.1};
    growFull = {false};
    fuzziness = {0, 0.1};
    steepness = {2};
    regularization = {0, 1};
  end
  
  methods (TestClassSetup)
    function setupClass(testCase)
      testCase.drawEnabled = false;
    end
  end
  
  methods (Test)
    function test0(testCase, fNum, m, ...
        minLeafSize, minGain, growFull, fuzziness, steepness, regularization)
      params = struct;
      params.modelSpec = '';
      params.minLeafSize = minLeafSize;
      params.minGain = minGain;
      params.growFull = growFull;
      params.fuzziness = fuzziness;
      params.steepness = steepness;
      params.regularization = regularization;
      testCase.reset(params, sprintf('_%02d_%03d', fNum, m));
      
      splitGainOptions = struct;
      splitGainOptions.minSize = minLeafSize;
      
      splits = {};
      splitOptions = struct;
      splitOptions.soft = fuzziness ~= 0;
      splitOptions.lambda = steepness;
      splits{1} = AxisSplit(splitOptions);
      
      treeModelOptions = struct;
      treeModelOptions.minGain = minGain;
      treeModelOptions.splitGain = GradientSplitGain(splitGainOptions);
      treeModelOptions.splits = splits;
      treeModelOptions.growFull = growFull;
      treeModelOptions.fuzziness = fuzziness;
      treeModelOptions.regularization = 0;
      treeModelFunc = @() GradientTreeModel(treeModelOptions);
      
      modelOptions = struct;
      modelOptions.nTrees = 10;
      modelOptions.treeFunc = treeModelFunc;
      modelFunc = @() XGBoostModel(modelOptions);
      
      [model, train, test, time] = testCase.testCoco(modelFunc, fNum, m);
    end
  end
end