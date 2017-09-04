classdef SSESplitGainTest < SplitGainTest  
  
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
      splitGain = SSESplitGain(options);
      
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
      options.degree = degree;
      options.polyMethod = polyMethod;
      splitGain = SSESplitGain(options);
      
      testMethod = strcat('split', testMethod);
      testCase.(testMethod)(splitGain);
    end
  end
end