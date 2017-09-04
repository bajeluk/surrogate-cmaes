classdef PolynomialModelTest < WeakModelTest
  
  properties (TestParameter)
    testMethod = {'ConstantFunction', 'LinearFunction', 'QuadraticFunction', ...
      'DependentFeatures', 'FewPoints'};
    modelSpec = {'constant', 'linear', 'quadratic'};
  end
  
  methods (Test)
    function test(testCase, testMethod, ...
        modelSpec)
      params = struct;
      params.modelSpec = modelSpec;
      testCase.reset(params, testMethod);
      
      modelOptions = struct;
      modelOptions.modelSpec = modelSpec;
      modelFunc = @() PolynomialModel(modelOptions);
      
      testMethod = strcat('test', testMethod);
      testCase.(testMethod)(modelFunc);
    end
  end
end