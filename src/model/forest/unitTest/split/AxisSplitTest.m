classdef AxisSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    nQuantize = {10, 100, -1};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        nQuantize, pca)
      params = struct;
      params.nQuantize = int2str(nQuantize);
      params.pca = int2str(pca);
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = AxisSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        nQuantize)
      params = struct;
      params.nQuantize = int2str(nQuantize);
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      split = AxisSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end