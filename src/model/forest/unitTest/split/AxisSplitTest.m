classdef AxisSplitTest < SplitTest
  methods (Test)
    
    function testOneLineParallelToAxis(testCase)
      % one line parallel with x axis
      n = 1000;
      m = 10;
      X = [rand(n, 1) * m, zeros(n, 1)];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testOneLineParallelToAxisSampled(testCase)
      % one line parallel with y axis
      n = 1000;
      m = 10;
      X = [zeros(n, 1), rand(n, 1) * m];
      % we expect the split to be in half
      y = (X(:, 2) <= m/2)*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testOneLine45Degrees(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 10;
      X = rand(n, 1) * m;
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testOneLine45DegreesPca(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 10;
      X = rand(n, 1) * m;
      X = [X X];
      % we expect the split to be in half
      y = (X(:, 1) <= m/2)*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = true;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testTwoParallelLines45Degrees(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 10;
      X = rand(n, 1) * m;
      half = (1:n)' <= n/2;
      X = [X(half, :)+1, X(half, :); ...
           X(~half, :),  X(~half, :)+1];
      % we expect to split the two parallel lines 
      y = half*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testTwoParallelLines45DegreesPca(testCase)
      rng('default');
      % one line where X1 == X2
      n = 1000;
      m = 10;
      X = rand(n, 1) * m;
      half = (1:n)' <= n/2;
      X = [X(half, :)+1, X(half, :); ...
           X(~half, :),  X(~half, :)+1];
      % we expect to split the two parallel lines 
      y = half*1;
      
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'linear';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = true;
      split = AxisSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    %{
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
    %}
  end
end