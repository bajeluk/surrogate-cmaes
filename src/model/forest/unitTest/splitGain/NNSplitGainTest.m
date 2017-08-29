classdef NNSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      splitGain = NNSplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      options = struct;
      splitGain = NNSplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      options = struct;
      splitGain = NNSplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

