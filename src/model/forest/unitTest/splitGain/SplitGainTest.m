classdef SplitGainTest < matlab.unittest.TestCase
  
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
        figure('Name', stack(2).name);
        scatter3(X(splitIdx, 1), X(splitIdx, 2), y(splitIdx));
        hold on;
        scatter3(X(~splitIdx, 1), X(~splitIdx, 2), y(~splitIdx));
        hold off;
        legend('yLeftWanted', 'yRightWanted');
      end
    end
    
    function [X, n, d, m] = generateRandom(testCase, n, d, m)
      rng('default');
      % random points
      if nargin < 2
        n = 1000; % number of points
      end
      if nargin < 3
        d = 2; % dimension
      end
      if nargin < 3
        m = 1; % range [-m, m]
      end
      X = (2 * rand(n, d) - 1) * m;
    end
    
    function [X, y, splitIdx, n, d, m] = generate(testCase, f1, f2)
      [X, n, d, m] = testCase.generateRandom();
      splitIdx = X(:, 1) <= 0;
      y1 = f1(X);
      if nargin < 3
        y2 = f1(X);
      else
        y2 = f2(X);
      end
      y = y1 .* splitIdx + y2 .* ~splitIdx;
      testCase.draw(X, y, splitIdx);
    end
  end
end

