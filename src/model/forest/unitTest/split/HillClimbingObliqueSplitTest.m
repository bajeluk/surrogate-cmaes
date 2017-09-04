classdef HillClimbingObliqueSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    nQuantize = {10, 100};
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
      splitOptions.nRepeats = 10;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = HillClimbingObliqueSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        nQuantize)
      params = struct;
      params.nQuantize = int2str(nQuantize);
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 10;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      split = HillClimbingObliqueSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
    
    function testQuadraticFeatures(testCase, testMethod, ...
        nQuantize)
      params = struct;
      params.nQuantize = int2str(nQuantize);
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 10;
      splitOptions.nQuantize = nQuantize;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.polynomial = 'quadratic';
      split = HillClimbingObliqueSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end