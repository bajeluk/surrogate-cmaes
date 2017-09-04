classdef DESplitGainTest < SplitGainTest
  
  properties (TestParameter)
    testMethod = {'Constant', 'Linear', 'Quadratic'};
    degree = {'constant', 'linear', 'quadratic'};
    polyMethod = {'default', 'regress'};
  end
  
  methods (Test)
    function test(testCase, testMethod)
      params = struct;
      testCase.reset(params, testMethod);
      
      options = struct;
      options.minSize = 5;
      splitGain = DESplitGain(options);
      
      testMethod = strcat('split', testMethod);
      testCase.(testMethod)(splitGain);
    end
    
    function testPoly(testCase, testMethod, ...
        degree, polyMethod)
      params = struct;
      params.degree = degree;
      params.polyMethod = polyMethod;
      testCase.reset(params, testMethod);
      
      options = struct;
      options.minSize = 5;
      options.degree = degree;
      options.polyMethod = polyMethod;
      splitGain = DESplitGain(options);
      
      testMethod = strcat('split', testMethod);
      testCase.(testMethod)(splitGain);
    end
  end
end

