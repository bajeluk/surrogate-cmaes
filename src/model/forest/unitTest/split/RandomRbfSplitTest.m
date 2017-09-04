classdef RandomRbfSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    metric = {'euclidean'};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        metric, pca)
      params = struct;
      params.metric = metric;
      params.pca = int2str(pca);
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.metric = metric;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = RandomRbfSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        metric)
      params = struct;
      params.metric = metric;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.metric = metric;
      splitOptions.transformationOptions = struct;
      split = RandomRbfSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end