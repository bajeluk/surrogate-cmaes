classdef RandomPolynomialTreeSplitGeneratorTest < matlab.unittest.TestCase
  methods (Test)
    function testOneLine(testCase)
      rng('default');
      nRepeats = 40;
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      % one variable should suffice
      generator = RandomPolynomialTreeSplitGenerator('linear', 1, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % it splits it aproximately in half
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testOneLine2(testCase)
      rng('default');
      nRepeats = 40;
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      % one variable should suffice
      generator = RandomPolynomialTreeSplitGenerator('linear2', 1, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % it splits it aproximately in half
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLines(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLines2(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear2', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % good split
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLinesCloser(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % not so good split
      verifyGreaterThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoParallelLinesCloser2(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear2', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % method 2 has good split already
      verifyLessThanOrEqual(testCase, best.err, n/100);
    end
    
    function testTwoParallelLinesCloserMoreRepeats(testCase)
      rng('default');
      nRepeats = 400;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % eventually good split
      verifyLessThanOrEqual(testCase, best.err, n/1000);
    end
    
    function testTwoCircles(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % quite good split
      verifyLessThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoCircles2(testCase)
      rng('default');
      nRepeats = 40;
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
      
      generator = RandomPolynomialTreeSplitGenerator('linear2', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % quite good split
      verifyLessThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoCirclesInside(testCase)
      rng('default');
      nRepeats = 40;
      % two circles
      n = 1000;
      m = 100;
      X = rand(n, 2) * 2 - 1;
      for i = 1:n
        while sqrt(sum(X(i, :).^2)) > 1
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = 1:n/2
        while sqrt(sum(X(i, :).^2)) < 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = n/2+1:n
        while sqrt(sum(X(i, :).^2)) > 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      half = (1:n)' <= n/2;
      % we expect to split the circles where they meet 
      y = half*1;
      
      generator = RandomPolynomialTreeSplitGenerator('linear', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % linear can't split it properly
      verifyGreaterThanOrEqual(testCase, best.err, n*3/10);
    end
    
    function testTwoCirclesInsideQuadratic(testCase)
      rng('default');
      nRepeats = 1000;
      % two circles
      n = 1000;
      m = 100;
      X = rand(n, 2) * 2 - 1;
      for i = 1:n
        while sqrt(sum(X(i, :).^2)) > 1
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = 1:n/2
        while sqrt(sum(X(i, :).^2)) < 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = n/2+1:n
        while sqrt(sum(X(i, :).^2)) > 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      half = (1:n)' <= n/2;
      % we expect to split the circles where they meet 
      y = half*1;
      
      generator = RandomPolynomialTreeSplitGenerator('quadratic', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % quadratic can split it quite ok
      verifyLessThanOrEqual(testCase, best.err, n*2/10);
    end
    
    function testTwoCirclesInsidePureQuadratic(testCase)
      rng('default');
      nRepeats = 1000;
      % two circles
      n = 1000;
      m = 100;
      X = rand(n, 2) * 2 - 1;
      for i = 1:n
        while sqrt(sum(X(i, :).^2)) > 1
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = 1:n/2
        while sqrt(sum(X(i, :).^2)) < 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      for i = n/2+1:n
        while sqrt(sum(X(i, :).^2)) > 0.5
          X(i, :) = rand(1, 2) * 2 - 1;
        end
      end
      half = (1:n)' <= n/2;
      % we expect to split the circles where they meet 
      y = half*1;
      
      generator = RandomPolynomialTreeSplitGenerator('purequadratic', 2, nRepeats);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      verifyEqual(testCase, count, nRepeats);
      % using only purequadratic decreases search space
      verifyLessThanOrEqual(testCase, best.err, n*2/10);
    end
  end
end