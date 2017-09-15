classdef RandomPolynomialSplitTest < SplitTest
  
  properties (TestParameter)
    testMethod = {'Flat', 'Axis', 'Linear', 'Linear2', 'Polynomial', ...
      'Circle', 'Atan', 'Parabola', 'Parabola2'};
    degree = {'linear', 'linear2', 'quadratic'};
    degreeQuadraticFeatures = {'linear', 'linear2'};
    pca = {false, true};
  end
  
  methods (Test)
    function testTwoLines(testCase, ...
        degree, pca)
      params = struct;
      params.degree = degree;
      params.pca = pca;
      testCase.reset(params);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.degree = degree;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = pca;
      split = RandomPolynomialSplit(splitOptions);
      
      [best] = testCase.splitTwoLines(split);
    end
    
    function test(testCase, testMethod, ...
        degree)
      params = struct;
      params.degree = degree;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.degree = degree;
      splitOptions.transformationOptions = struct;
      split = RandomPolynomialSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
    
    function testQuadraticFeatures(testCase, testMethod, ...
        degreeQuadraticFeatures)
      params = struct;
      params.degree = degreeQuadraticFeatures;
      testCase.reset(params, testMethod);
      
      splitOptions = struct;
      splitOptions.nRepeats = 1000;
      splitOptions.degree = degreeQuadraticFeatures;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.polynomial = 'quadratic';
      split = RandomPolynomialSplit(splitOptions);
      
      testMethod = strcat('split', testMethod);
      [best] = testCase.(testMethod)(split);
    end
  end
end