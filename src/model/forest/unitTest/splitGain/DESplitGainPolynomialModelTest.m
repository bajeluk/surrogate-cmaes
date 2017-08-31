classdef DESplitGainPolynomialModelTest < SplitGainTest
  
  methods (Test)
   function testConstantConstantPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'constant';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearConstantPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'constant';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticConstantPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'constant';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantLinearPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'linear';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearLinearPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'linear';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticLinearPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'linear';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantQuadraticPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'quadratic';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearQuadraticPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'quadratic';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticQuadraticPolynomialModel(testCase)
      modelOptions = struct;
      modelOptions.modelSpec = 'quadratic';
      options = struct;
      options.modelFunc = @() PolynomialModel(modelOptions);
      splitGain = DESplitGain(options);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

