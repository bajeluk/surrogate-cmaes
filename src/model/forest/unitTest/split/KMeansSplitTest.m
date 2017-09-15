classdef KMeansSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    discrimType = {'linear', 'quadratic'};
    includeInput = {false, true};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        discrimType, includeInput, pca)
      params = struct;
      params.discrimType = discrimType;
      params.includeInput = includeInput;
      params.pca = pca;
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nRepeats = 10;
      splitOptions.discrimType = discrimType;
      splitOptions.includeInput = includeInput;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = KMeansSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        discrimType, includeInput)
      params = struct;
      params.discrimType = discrimType;
      params.includeInput = includeInput;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 10;
      splitOptions.discrimType = discrimType;
      splitOptions.includeInput = includeInput;
      splitOptions.transformationOptions = struct;
      split = KMeansSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end