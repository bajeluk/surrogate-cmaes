classdef ConstantModelTest < WeakModelTest
  
  properties (TestParameter)
    testMethod = {'ConstantFunction', 'LinearFunction', 'QuadraticFunction', ...
      'DependentFeatures', 'FewPoints'};
    modelSpec = {'constant', 'linear', 'quadratic'};
  end
  
  methods (Test)
    function testConstantFunctionSet(testCase)
      params = struct;
      testCase.reset(params);
      
      modelOptions = struct('coeff', 0.75);
      modelFunc = @() ConstantModel(modelOptions);
      
      testCase.testConstantFunction(modelFunc);
    end
    
    function test(testCase, testMethod)
      params = struct;
      testCase.reset(params, testMethod);
      
      modelOptions = struct;
      modelFunc = @() ConstantModel(modelOptions);
      
      testMethod = strcat('test', testMethod);
      testCase.(testMethod)(modelFunc);
    end
  end
end