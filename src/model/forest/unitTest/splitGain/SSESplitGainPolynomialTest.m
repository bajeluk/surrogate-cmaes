classdef SSESplitGainPolynomialTest < SplitGainTest
  
  methods (Test)
    function testConstantConstantModel(testCase)
      options = struct;
      options.degree = 'constant';
      splitGain = SSESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearConstantModel(testCase)
      options = struct;
      options.degree = 'constant';
      splitGain = SSESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticConstantModel(testCase)
      options = struct;
      options.degree = 'constant';
      splitGain = SSESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantLinearModel(testCase)
      options = struct;
      options.degree = 'linear';
      splitGain = SSESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearLinearModel(testCase)
      options = struct;
      options.degree = 'linear';
      splitGain = SSESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticLinearModel(testCase)
      options = struct;
      options.degree = 'linear';
      splitGain = SSESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantQuadraticModel(testCase)
      options = struct;
      options.degree = 'quadratic';
      splitGain = SSESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearQuadraticModel(testCase)
      options = struct;
      options.degree = 'quadratic';
      splitGain = SSESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticQuadraticModel(testCase)
      options = struct;
      options.degree = 'quadratic';
      splitGain = SSESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

