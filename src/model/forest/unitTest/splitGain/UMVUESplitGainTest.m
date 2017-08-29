classdef UMVUESplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      splitGain = UMVUESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      options = struct;
      splitGain = UMVUESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      options = struct;
      splitGain = UMVUESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

