classdef GMMTreeSplitGeneratorTest < matlab.unittest.TestCase
  methods (Test)
    function testOneGaussian(testCase)
      rng('default');
      nRepeats = 2;
      % one gaussian in one line
      n = 1000;
      X = randn(n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= 0)*1;
      
      generator = GMMTreeSplitGenerator(nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % can't really split one gaussian
      verifyGreaterThanOrEqual(testCase, best.err, n/10);
    end
    
    function testTwoGaussians(testCase)
      rng('default');
      nRepeats = 2;
      % two gaussians in one line
      n = 1000;
      m = 10;
      X = [[randn(n/2, 1); randn(n/2, 1) + m]];
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = GMMTreeSplitGenerator(nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelGaussians(testCase)
      rng('default');
      nRepeats = 2;
      % two gaussians in two parallel lines
      n = 1000;
      m = 10;
      X = randn(n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - m/2;
      X(~half, 2) = X(~half, 2) + m/2;
      % we expect to split the two parallel gaussians 
      y = half*1;
      
      generator = GMMTreeSplitGenerator(nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
  end
end