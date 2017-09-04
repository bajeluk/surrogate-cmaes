classdef GradientSplitGainTest < SplitGainTest
  
  properties (TestParameter)
    testMethod = {'Constant', 'Linear', 'Quadratic'};
    regularization = {0, 1};
  end
  
  methods (Test)
    function test(testCase, testMethod, ...
        regularization)
      params = struct;
      params.regularization = int2str(regularization);
      testCase.reset(params, testMethod);
      
      options = struct;
      options.weightedGain = false;
      options.regularization = regularization;
      splitGain = GradientSplitGain(options);
      g1 = @(y, yPred) 2*(yPred - y);
      g2 = @(y, yPred) y*0 + 2;
      g = @(y) [g1(y, 0), g2(y, 0)];
      
      testMethod = strcat('split', testMethod);
      testCase.(testMethod)(splitGain, g);
    end
  end
end

