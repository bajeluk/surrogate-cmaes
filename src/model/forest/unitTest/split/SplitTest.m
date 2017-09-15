classdef (Abstract) SplitTest < Test
  
  methods (Access = protected)
    function [X, y] = generate(testCase, f, phi)
      d = 2;
      n = 250 * d;
      m = 100;
      X = testCase.generateInput(n, d, 0, m);
      mu = [m/2 m/2];
      if nargin < 3
        phi = pi/4;
      end
      R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
      XT = bsxfun(@minus, X, mu) * R;
      y = f(XT, m);
    end
    
    function [best] = findBest(testCase, X, y, split, splitGain)
      splitGain = splitGain.reset(X, y);
      split = split.reset(X, y);
      best = split.get(splitGain);
      idx = best.splitter(X) <= 0.5;
      testCase.drawSplit(X, y, idx);
    end
    
    function [best] = splitTwoLines(testCase, split, splitGain)
      % one line where X1 == X2
      d = 2;
      n = 250 * d;
      m = 10;
      X = rand(n, 1) * m;
      half = (1:n)' <= n/2;
      X = [X(half, :)+1, X(half, :); ...
           X(~half, :),  X(~half, :)+1];
      % we expect to split the two parallel lines 
      y = half*1;
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitFlat(testCase, split, splitGain)
      f = @(X, m) X(:,1) * 0;
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitAxis(testCase, split, splitGain)
      f = @(X, m) X(:,1) > 0;
      [X, y] = testCase.generate(f, 0);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitLinear(testCase, split, splitGain)
      f = @(X, m) X(:,1) > 0;
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitLinear2(testCase, split, splitGain)
      f = @(X, m) (X(:,1) > -m/4) .* (X(:,1) < m/4);
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitPolynomial(testCase, split, splitGain)
      f = @(X, m) X(:,1).^2 > X(:,2) * m;
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitCircle(testCase, split, splitGain)
      f = @(X, m) X(:,1).^2 + X(:,2).^2 < (m/4)^2;
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitAtan(testCase, split, splitGain)
      f = @(X, m) atan(X(:,1));
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitParabola(testCase, split, splitGain)
      f = @(X, m) (-X(:,1).^2 + (m/4)^2) .* (X(:,1).^2 < (m/4)^2);
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGainOptions.degree = 'quadratic';
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitParabola2(testCase, split, splitGain)
      f = @(X, m) (-(X(:,1).^2 + X(:,2).^2) + (m/4)^2) .* (X(:,1).^2 + X(:,2).^2 < (m/4)^2);
      [X, y] = testCase.generate(f);
      if nargin < 3
        splitGainOptions = struct;
        splitGainOptions.minSize = 5;
        splitGainOptions.degree = 'quadratic';
        splitGain = SSESplitGain(splitGainOptions);
      end
      best = testCase.findBest(X, y, split, splitGain);
    end
  end
end

