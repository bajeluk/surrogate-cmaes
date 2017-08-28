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
    
    function testConstant10(testCase)
      splitGain = NNSplitGain(10);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear10(testCase)
      splitGain = NNSplitGain(10);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic10(testCase)
      splitGain = NNSplitGain(10);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

