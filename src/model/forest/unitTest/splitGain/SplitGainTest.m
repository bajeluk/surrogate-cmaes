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
        figure('Name', stack(4).name);
        scatter3(X(splitIdx, 1), X(splitIdx, 2), y(splitIdx));
        hold on;
        scatter3(X(~splitIdx, 1), X(~splitIdx, 2), y(~splitIdx));
        hold off;
        legend('yLeftWanted', 'yRightWanted');
      end
    end
    
    function [X, n, d, minVal, maxVal] = ...
        generateRandom(testCase, n, d, minVal, maxVal)
      rng('default');
      % random points
      if nargin < 2
        n = 1000; % number of points
      end
      if nargin < 3
        d = 2; % dimension
      end
      if nargin < 3
        minVal = -100; % range [minVal, maxVal]
      end
      if nargin < 4
        maxVal = 100; % range [minVal, maxVal]
      end
      X = minVal + (maxVal - minVal) * rand(n, d);
    end
    
    function [X, y, splitIdx, n, d, minVal, maxVal, splitVal] = ...
        generate(testCase, f1, f2, splitPercentage)
      [X, n, d, minVal, maxVal] = testCase.generateRandom();
      if nargin < 4
        splitPercentage = 0.5;
      end
      splitVal = minVal * (1-splitPercentage) + maxVal * splitPercentage;
      splitIdx = X(:, 1) < splitVal;
      y1 = f1(X);
      if nargin < 3
        y2 = f1(X);
      else
        y2 = f2(X);
      end
      y = y1 .* splitIdx + y2 .* ~splitIdx;
      testCase.draw(X, y, splitIdx);
    end
    
    function [values, maxGain, maxGainSplit] = ...
        testGains(testCase, splitGain, f1, f2, splitPercentage)
      if nargin < 5
        splitPercentage = 0.5;
      end
      [X, y, splitIdx, n, d, minVal, maxVal, midVal] = ...
        testCase.generate(f1, f2, splitPercentage);
      splitGain = splitGain.reset(X, y);
      splitPoints = unique(X(:, 1));
      values = -inf(numel(splitPoints), 1);
      for i = 1:numel(splitPoints)
        x = splitPoints(i);
        splitter = @(X) X(:, 1) < x;
        values(i) = splitGain.get(splitter);
      end
      [maxGain, maxGainPos] = max(values);
      maxGainSplit = splitPoints(maxGainPos);
      
      if testCase.drawEnabled
        stack = dbstack;
        figure('Name', stack(2).name);
        plot(splitPoints, values);
      end
    end
    
    function testAxisConstant(testCase, splitGain)
      f = @(X) repmat(max(X(:, 1)), size(X, 1), 1) + randn(size(X, 1), 1) * 0.001;
      testCase.testAxis(splitGain, f);
    end
    
    function testAxisLinear(testCase, splitGain)
      f = @(X) X(:, 1) + randn(size(X, 1), 1) * 0.001;
      testCase.testAxis(splitGain, f);
    end
    
    function testAxisQuadratic(testCase, splitGain)
      f = @(X) X(:, 1).^2 + randn(size(X, 1), 1) * 0.001;
      testCase.testAxis(splitGain, f);
    end
    
    function testAxis(testCase, splitGain, f)
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, f);
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, @(X)-f(X));
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, @(X)-f(X), 0.75);
    end
  end
end

