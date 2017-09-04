classdef (Abstract) SplitGainTest < Test
  
  methods (Access = protected)    
    function [X, y, splitIdx, minVal, maxVal, splitVal] = ...
        generate(testCase, f1, f2, splitPercentage)
      [X, minVal, maxVal] = testCase.generateInput();
      if nargin < 4
        splitPercentage = 0.5;
      end
      splitVal = minVal * (1-splitPercentage) + maxVal * splitPercentage;
      splitIdx = X(:, 1) < splitVal;
      y1 = f1(X - splitVal);
      if nargin < 3
        y2 = f1(X - splitVal);
      else
        y2 = f2(X - splitVal);
      end
      y = y1 .* splitIdx + y2 .* ~splitIdx;
      %y = y + randn(size(y));
    end
    
    function [values, maxGain, maxGainSplit] = ...
        findBest(testCase, splitGain, f1, f2, g, splitPercentage)
      rng('default');
      if nargin < 6
        splitPercentage = 0.5;
      end
      [X, y, splitIdx, minVal, maxVal, midVal] = ...
        testCase.generate(f1, f2, splitPercentage);
      splitGain = splitGain.reset(X, g(y));
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
        splitIdx = X(:, 1) < splitPoints(maxGainPos);
        testCase.iPlot = testCase.iPlot + 1;
        subplot(3, 2, testCase.iPlot);
        testCase.drawSplit(X, y, splitIdx);
        title('split');
        testCase.iPlot = testCase.iPlot + 1;
        subplot(3, 2, testCase.iPlot);
        plot(splitPoints, values);
        title('gain');
      end
    end
    
    function splitConstant(testCase, splitGain, g)
      if nargin < 3
        g = @(y) y;
      end
      f = @(X) repmat(100, size(X, 1), 1) + randn(size(X, 1), 1) * 0.01;
      testCase.split(splitGain, f, g);
    end
    
    function splitLinear(testCase, splitGain, g)
      if nargin < 3
        g = @(y) y;
      end
      f = @(X) X(:, 1) + randn(size(X, 1), 1) * 0.01;
      testCase.split(splitGain, f, g);
    end
    
    function splitQuadratic(testCase, splitGain, g)
      if nargin < 3
        g = @(y) y;
      end
      f = @(X) X(:, 1).^2 + randn(size(X, 1), 1) * 0.01;
      testCase.split(splitGain, f, g);
      %{
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, f, g);
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, f, g, 0.75);
      [values, maxGain, maxGainSplit] = ...
        testCase.testGains(splitGain, f, @(X)-f(X), g);
      %}
    end
    
    function split(testCase, splitGain, f, g)
      [values, maxGain, maxGainSplit] = ...
        testCase.findBest(splitGain, f, f, g);
      [values, maxGain, maxGainSplit] = ...
        testCase.findBest(splitGain, f, @(X)-f(X), g);
      [values, maxGain, maxGainSplit] = ...
        testCase.findBest(splitGain, f, @(X)-f(X), g, 0.75);
    end
  end
end

