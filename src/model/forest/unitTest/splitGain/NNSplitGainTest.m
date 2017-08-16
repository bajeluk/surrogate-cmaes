classdef NNSplitGainTest < SplitGainTest
  
  methods (Test)
    function testConstant1(testCase)
      f = @(X) ones(size(X, 1), 1);
      [X, y, splitIdx, n, d, m] = generate(testCase, f, f);
      
      splitGain = NNSplitGain();
      splitGain = splitGain.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/2;
      gain = splitGain.eval(splitter);
      % no gain
      testCase.verifyLessThanOrEqual(gain, -realmax/2);
      
      % try to split in half
      splitter = @(X) X(:, 1) <= 0;
      gain = splitGain.eval(splitter);
      % no gain
      testCase.verifyLessThanOrEqual(gain, -realmax/2);
    end
    
    function testConstant2(testCase)
      f = @(X) ones(size(X, 1), 1);
      [X, y, splitIdx, n, d, m] = generate(testCase, f, @(X)-f(X));
      
      splitGain = NNSplitGain();
      splitGain = splitGain.reset(X, y);
      
      % try to split in quarter
      splitter = @(X) X(:, 1) <= m/2;
      gain = splitGain.eval(splitter);
      % no gain
      testCase.verifyLessThanOrEqual(gain, -realmax/2);
      
      % try to split in half
      splitter = @(X) X(:, 1) <= 0;
      gain = splitGain.eval(splitter);
      % big gain
      testCase.verifyGreaterThanOrEqual(gain, realmax-1);
    end
  end
end

