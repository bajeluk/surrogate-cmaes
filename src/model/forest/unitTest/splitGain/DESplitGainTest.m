classdef DESplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      splitGain = DESplitGain();
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinear(testCase)
      splitGain = DESplitGain();
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadratic(testCase)
      splitGain = DESplitGain();
      testCase.testAxisQuadratic(splitGain);
    end
    
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
    
    function testConstantConstantPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'constant');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearConstantPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'constant');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticConstantPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'constant');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantLinearPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearLinearPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticLinearPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'linear');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisQuadratic(splitGain);
    end
    
    function testConstantQuadraticPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'quadratic');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisConstant(splitGain);
    end
    
    function testLinearQuadraticPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'quadratic');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisLinear(splitGain);
    end
    
    function testQuadraticQuadraticPolynomialModel(testCase)
      modelOptions = struct('modelSpec', 'quadratic');
      modelFunc = @(xMean) PolynomialModel(modelOptions, xMean);
      splitGain = DESplitGain(modelFunc);
      testCase.testAxisQuadratic(splitGain);
    end
  end
end

