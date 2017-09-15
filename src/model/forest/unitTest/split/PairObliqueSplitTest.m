classdef PairObliqueSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    nQuantize = {10, -1};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        nQuantize, pca)
      params = struct;
      params.nQuantize = nQuantize;
      params.pca = pca;
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.nValues = floor(sqrt(1000));
      splitOptions.transformationOptions.pca = pca;
      split = PairObliqueSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        nQuantize)
      params = struct;
      params.nQuantize = nQuantize;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.nValues = floor(sqrt(1000));
      split = PairObliqueSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end