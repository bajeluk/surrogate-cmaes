classdef RandomSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Random', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
      pca)
      params = struct;
      params.pca = pca;
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = RandomSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod)
      params = struct;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.transformationOptions = struct;
      split = RandomSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end