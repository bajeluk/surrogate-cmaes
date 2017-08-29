classdef VarSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      splitGain = VarSplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      options = struct;
      splitGain = VarSplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      options = struct;
      splitGain = VarSplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

