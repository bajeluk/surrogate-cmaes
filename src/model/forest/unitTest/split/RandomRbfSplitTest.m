classdef RandomRbfTreeSplitGeneratorTest < matlab.unittest.TestCase
  methods (Test)
    function testOneLine(testCase)
      rng('default');
      nRepeats = 100;
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = RandomRbfTreeSplitGenerator('euclidean', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % it splits it aproximately in half eventually
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoLines(testCase)
      rng('default');
      nRepeats = 100;
      % two gaussians in one line
      n = 1000;
      m = 100;
      X = [[randi(m, n/2, 1); randi(m, n/2, 1) + 2*m]];
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= 1.5*m)*1;
      
      generator = RandomRbfTreeSplitGenerator('euclidean', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/1000);
    end
    
    function testTwoParallelLines(testCase)
      rng('default');
      nRepeats = 100;
      % two parallel lines
      n = 1000;
      m = 100;
      X = randi(100, n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - 60;
      X(~half, 2) = X(~half, 2) + 60;
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = RandomRbfTreeSplitGenerator('euclidean', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLinesCloser(testCase)
      rng('default');
      nRepeats = 100;
      % two parallel lines
      n = 1000;
      m = 100;
      X = randi(100, n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - 10;
      X(~half, 2) = X(~half, 2) + 10;
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = RandomRbfTreeSplitGenerator('euclidean', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % nt so good split
      verifyGreaterThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoParallelLinesCloserMahalanobis(testCase)
      rng('default');
      nRepeats = 100;
      % two parallel lines
      n = 1000;
      m = 100;
      X = randi(100, n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - 10;
      X(~half, 2) = X(~half, 2) + 10;
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = RandomRbfTreeSplitGenerator('mahalanobis', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % mahalanobis is better than euclidean
      verifyLessThanOrEqual(testCase, best.err, n*1/100);
    end
    
    function testTwoCircles(testCase)
      rng('default');
      nRepeats = 100;
      % two circles
      n = 1000;
      m = 100;
      X = rand(n, 2) * 2 - 1;
      for i = 1:n
        while sqrt(sum(X(i, :).^2)) > 1
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      half = (1:n)' <= n/2;
      X(~half, :) = X(~half, :) + 1;
      % we expect to split the circles where they meet
      y = half*1;
      
      generator = RandomRbfTreeSplitGenerator('euclidean', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % quite good split
      verifyLessThanOrEqual(testCase, best.err, n/10);
    end
  end
end