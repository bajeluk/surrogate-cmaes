classdef DESplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      splitGain = DESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      options = struct;
      splitGain = DESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      options = struct;
      splitGain = DESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

