classdef DEMSDSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      splitGain = DEMSDSplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      options = struct;
      splitGain = DEMSDSplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      options = struct;
      splitGain = DEMSDSplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

