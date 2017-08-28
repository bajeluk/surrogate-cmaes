classdef NNSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      splitGain = NNSplitGain();
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      splitGain = NNSplitGain();
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      splitGain = NNSplitGain();
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

