classdef DESplitGainPolynomialTest < SplitGainTest
  
  methods (Test)
    function testConstantConstantModel(testCase)
      splitGain = DESplitGain('constant');
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearConstantModel(testCase)
      splitGain = DESplitGain('constant');
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticConstantModel(testCase)
      splitGain = DESplitGain('constant');
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantLinearModel(testCase)
      splitGain = DESplitGain('linear');
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearLinearModel(testCase)
      splitGain = DESplitGain('linear');
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticLinearModel(testCase)
      splitGain = DESplitGain('linear');
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantQuadraticModel(testCase)
      splitGain = DESplitGain('quadratic');
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearQuadraticModel(testCase)
      splitGain = DESplitGain('quadratic');
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticQuadraticModel(testCase)
      splitGain = DESplitGain('quadratic');
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

