classdef AxisTreeSplitGeneratorTest < matlab.unittest.TestCase
  methods (Test)
    function testOneLineParallelToAxis(testCase)
      rng('default');
      % one line parallel with x axis
      n = 1000;
      m = 100;
      X = [randi(m, n, 1), zeros(n, 1)];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = AxisTreeSplitGenerator(1, 1, false);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, m+1);
      % perfect split in the middle
      verifyEqual(testCase, best.err, 0);
      verifyEqual(testCase, vars.feature, 1);
      verifyEqual(testCase, vars.treshold, m/2);
    end
    
    function testOneLineParallelToAxisSampled(testCase)
      rng('default');
      % one line parallel with y axis
      n = 1000;
      m = 100;
      X = [zeros(n, 1), randi(m, n, 1)];
      % we expect the split to be in half
      y = (X(:, 2) <= m/2)*1;
      
      generator = AxisTreeSplitGenerator(1, 0.5, false);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, m/2+1);
      % perfect split in the middle
      verifyEqual(testCase, best.err, 0);
      verifyEqual(testCase, vars.feature, 2);
      verifyEqual(testCase, vars.treshold, m/2);
    end
    
    function testOneLine45Degrees(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = AxisTreeSplitGenerator(1, 1, false);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, 2*m);
      % perfect split in the middle
      verifyEqual(testCase, best.err, 0);
      verifyEqual(testCase, vars.treshold, m/2);
    end
    
    function testOneLine45DegreesPca(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      generator = AxisTreeSplitGenerator(1, 1, true);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, m+1);
      % perfect split in the middle
      verifyEqual(testCase, best.err, 0);
      % PCA rotates the line so that it's parallel with X1 axis
      verifyEqual(testCase, vars.feature, 1);
      verifyLessThan(testCase, abs(vars.treshold - 70.7107), 0.001);
    end
    
    function testTwoParallelLines45Degrees(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      half = (1:n)' <= n/2;
      X = [X(half, :)+1, X(half, :); ...
           X(~half, :),  X(~half, :)+1];
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = AxisTreeSplitGenerator(1, 1, false);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, 2*(m+1));
      % bad split
      verifyGreaterThan(testCase, best.err, n * 4/10);
    end
    
    function testTwoParallelLines45DegreesPca(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 100;
      X = randi(m, n, 1);
      half = (1:n)' <= n/2;
      X = [X(half, :)+1, X(half, :); ...
           X(~half, :),  X(~half, :)+1];
      % we expect to split the two parallel lines 
      y = half*1;
      
      generator = AxisTreeSplitGenerator(1, 1, true);
      [best, count, vars] = TreeSplitGeneratorTest.findBest(X, y, generator);
      
      % only distinct values are considered
      verifyLessThanOrEqual(testCase, count, 2*n);
      % bad split
      verifyEqual(testCase, best.err, 0);
      % PCA rotates the line so that it's parallel with X1 axis
      verifyEqual(testCase, vars.feature, 2);
      verifyLessThan(testCase, abs(vars.treshold - 0), 1);
    end
    
    function testTime(testCase)
      rng('default');
      X = randi(100, 1000, 20);
      y = randi(100, size(X, 1));

      tic
      n = 0;
      for d = 1:size(X, 2)
        for i = 1:size(X, 1)
          t = X(i, d);
          f = @(X) X(i, d) <= t;
          n = n+1;
        end
      end
      nLoops = n;
      tLoops = toc

      tic
      n = 0;
      for d = 1:size(X, 2)
        x = unique(X(:, d));
        for i = 1:size(x, 1)
          t = x(i, 1);
          f = @(X) X(i, d) <= t;
          n = n+1;
        end
      end
      nLoopsUnique = n;
      tLoopsUnique = toc

      tic
      n = 0;
      it = AxisTreeSplitGenerator(1, 1, false);
      it.reset(X, y);
      while it.hasNext()
        f = it.next();
        n = n+1;
      end
      nIterator = n;
      tIterator = toc

      tic
      n = 0;
      it = AxisTreeSplitGenerator(1, 1, true);
      it.reset(X, y);
      while it.hasNext()
        f = it.next();
        n = n+1;
      end
      nIteratorPca = n;
      tIteratorPca = toc

      % iterator is comparable with loops
      verifyEqual(testCase, nIterator, nLoopsUnique);
      verifyLessThan(testCase, abs(tIterator-tLoopsUnique), 0.1);
      
      verifyEqual(testCase, nIteratorPca, nLoops);
      verifyLessThan(testCase, abs(tIteratorPca-tLoops), 1);
    end
  end
end