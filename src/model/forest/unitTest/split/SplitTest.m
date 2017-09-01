classdef SplitTest < matlab.unittest.TestCase
  
  properties
    drawEnabled = true
  end
  
  methods (TestMethodSetup)
    function setup(testCase)
      rng('default');
    end
  end

  methods (TestMethodTeardown)
    function teardown(testCase)
    end
  end
  
  methods (Access = protected)
    function draw(testCase, X, y, splitIdx)      
      if size(X, 2) == 2 && testCase.drawEnabled
        stack = dbstack;
        name = strcat(stack(3).name);
        figure('Name', name);
        %scatter(X(splitIdx, 1), y(splitIdx));
        scatter3(X(splitIdx, 1), X(splitIdx, 2), y(splitIdx));
        hold on;
        %scatter(X(~splitIdx, 1), y(~splitIdx));
        scatter3(X(~splitIdx, 1), X(~splitIdx, 2), y(~splitIdx));
        hold off;
        legend('left', 'right');
      end
    end
    
    function [best] = findBest(testCase, X, y, split, splitGain)
      splitGain = splitGain.reset(X, y);
      split = split.reset(X, y);
      best = split.get(splitGain);
      idx = best.splitter(X);
      testCase.draw(X, y, idx);
    end
    
    function [X, y] = generate(testCase, f, phi)
      n = 1000;
      m = 100;
      X = rand(n, 2) * m;
      mu = [m/2 m/2];
      if nargin < 3
        phi = pi/4;
      end
      R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
      XT = (X - mu) * R;
      y = f(XT, m);
    end
    
    function [best] = splitFlat(testCase, split, splitGain)
      f = @(X, m) X(:,1) * 0;
      [X, y] = testCase.generate(f);
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitAxis(testCase, split, splitGain)
      f = @(X, m) X(:,1) > 0;
      [X, y] = testCase.generate(f, 0);
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitLine(testCase, split, splitGain)
      f = @(X, m) X(:,1) > 0;
      [X, y] = testCase.generate(f);
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitTwoLines(testCase, split, splitGain)
      f = @(X, m) (X(:,1) > -m/4) .* (X(:,1) < m/4);
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitPolynomial(testCase, split, splitGain)
      f = @(X, m) X(:,1).^2 > X(:,2) * m;
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitCircle(testCase, split, splitGain)
      f = @(X, m) X(:,1).^2 + X(:,2).^2 < (m/4)^2;
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitAtan(testCase, split, splitGain)
      f = @(X, m) atan(X(:,1));
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitParabola(testCase, split, splitGain)
      f = @(X, m) (-X(:,1).^2 + (m/4)^2) .* (X(:,1).^2 < (m/4)^2);
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
    
    function [best] = splitParabola2(testCase, split, splitGain)
      f = @(X, m) (-(X(:,1).^2 + X(:,2).^2) + (m/4)^2) .* (X(:,1).^2 + X(:,2).^2 < (m/4)^2);
      [X, y] = testCase.generate(f);      
      best = testCase.findBest(X, y, split, splitGain);
    end
  end
end

