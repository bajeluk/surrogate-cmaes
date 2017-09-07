classdef ResidualObliqueSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    modelSpec = {'constant', 'linear', 'quadratic'};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        modelSpec, pca)
      params = struct;
      params.modelSpec = modelSpec;
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.modelSpec = modelSpec;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = ResidualObliqueSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        modelSpec)
      params = struct;
      params.modelSpec = modelSpec;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.modelSpec = modelSpec;
      splitOptions.transformationOptions = struct;
      split = ResidualObliqueSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end