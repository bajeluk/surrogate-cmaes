classdef RandomRbfSplitChebychevTest < SplitTest
  methods (Test)
    
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
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
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
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.transformationOptions.pca = true;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.findBest(X, y, split, splitGain);
    end
    
    function testFlat(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitFlat(split, splitGain);
    end
    
    function testAxis(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitAxis(split, splitGain);
    end
    
    function testLine(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitLine(split, splitGain);
    end
    
    function testTwoLines(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitTwoLines(split, splitGain);
    end
    
    function testPolynomial(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitPolynomial(split, splitGain);
    end
    
    function testCircle(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'constant';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitCircle(split, splitGain);
    end
    
    function testAtan(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'quadratic';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitAtan(split, splitGain);
    end
    
    function testParabola(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'quadratic';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitParabola(split, splitGain);
    end
    
    function testParabola2(testCase)
      splitGainOptions = struct;
      splitGainOptions.minSize = 5;
      splitGainOptions.degree = 'quadratic';
      splitGain = SSESplitGain(splitGainOptions);
      splitOptions = struct;
      splitOptions.transformationOptions = struct;
      splitOptions.nRepeats = 100;
      splitOptions.metric = 'chebychev';
      split = RandomRbfSplit(splitOptions);
      [best] = testCase.splitParabola2(split, splitGain);
    end
  end
end