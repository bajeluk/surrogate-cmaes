classdef UMVUESplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      splitGain = UMVUESplitGain();
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      splitGain = UMVUESplitGain();
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      splitGain = UMVUESplitGain();
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

