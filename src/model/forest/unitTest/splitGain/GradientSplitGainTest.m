classdef GradientSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant(testCase)
      options = struct;
      options.weightedGain = false;
      options.regularization = 0;
      splitGain = GradientSplitGain(options);
      g1 = @(y, yPred) 2*(yPred - y);
      g2 = @(y, yPred) y*0 + 2;
      g = @(y) [g1(y, 0), g2(y, 0)];
      testCase.testAxisConstant(splitGain, g);
    end
    
    function testLinear(testCase)
      options = struct;
      options.weightedGain = false;
      options.regularization = 0;
      splitGain = GradientSplitGain(options);
      g1 = @(y, yPred) 2*(yPred - y);
      g2 = @(y, yPred) y*0 + 2;
      g = @(y) [g1(y, 0), g2(y, 0)];
      testCase.testAxisLinear(splitGain, g);
    end
    
    function testQuadratic(testCase)
      options = struct;
      options.weightedGain = false;
      options.regularization = 0;
      splitGain = GradientSplitGain(options);
      g1 = @(y, yPred) 2*(yPred - y);
      g2 = @(y, yPred) y*0 + 2;
      g = @(y) [g1(y, 0), g2(y, 0)];
      testCase.testAxisQuadratic(splitGain, g);
    end
  end
end

