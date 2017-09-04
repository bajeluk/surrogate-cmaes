classdef DENNSplitGainTest < SplitGainTest

  properties (TestParameter)
    testMethod = {'Constant', 'Linear', 'Quadratic'};
  end
  
  methods (Test)
    function test(testCase, testMethod)
      params = struct;
      testCase.reset(params, testMethod);
      
      options = struct;
      splitGain = DENNSplitGain(options);

      testMethod = strcat('split', testMethod);
      testCase.(testMethod)(splitGain);
    end
  end
end

