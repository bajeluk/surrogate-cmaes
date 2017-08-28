classdef VarSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      splitGain = VarSplitGain();
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      splitGain = VarSplitGain();
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      splitGain = VarSplitGain();
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

