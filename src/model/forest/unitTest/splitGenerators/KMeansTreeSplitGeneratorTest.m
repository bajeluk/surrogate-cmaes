classdef KMeansTreeSplitGeneratorTest < matlab.unittest.TestCase
  methods (Test)
    function testOneLine(testCase)
      rng('default');
      nRepeats = 2;
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', false, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % it splits it aproximately in half
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoLines(testCase)
      rng('default');
      nRepeats = 2;
      % two gaussians in one line
      n = 1000;
      m = 100;
      X = [[randi(m, n/2, 1); randi(m, n/2, 1) + 2*m]];
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= 1.5*m)*1;
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', false, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/1000);
    end
    
    function testTwoParallelLines(testCase)
      rng('default');
      nRepeats = 2;
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
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', false, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLinesCloser(testCase)
      rng('default');
      nRepeats = 2;
      % two parallel lines
      n = 1000;
      m = 100;
      X = randi(100, n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - 40;
      X(~half, 2) = X(~half, 2) + 40;
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', false, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % not so good split
      verifyGreaterThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoParallelLinesCloserNormalized(testCase)
      rng('default');
      nRepeats = 2;
      % two parallel lines
      n = 1000;
      m = 100;
      X = randi(100, n, 1);
      X = [X X];
      half = (1:n)' <= n/2;
      X(half, 2) = X(half, 2) - 40;
      X(~half, 2) = X(~half, 2) + 40;
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', true, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % even worse split
      verifyGreaterThanOrEqual(testCase, best.err, n*4/10);
    end
    
    function testTwoCircles(testCase)
      rng('default');
      nRepeats = 2;
      % two parallel lines
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
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = KMeansTreeSplitGenerator('sqeuclidean', false, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % quite good split
      verifyLessThanOrEqual(testCase, best.err, n/10);
    end
  end
end