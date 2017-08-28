classdef DESplitGainPolynomialModelTest < SplitGainTest
  
  methods (Test)
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

