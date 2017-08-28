classdef DifferentialEntropySplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      splitGain = DifferentialEntropySplitGain();
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      splitGain = DifferentialEntropySplitGain();
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      splitGain = DifferentialEntropySplitGain();
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

