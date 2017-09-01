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
    
    function [best] = ...
        findBest(testCase, X, y, split, splitGain)
      splitGain = splitGain.reset(X, y);
      split = split.reset(X, y);
      best = split.get(splitGain);
      idx = best.splitter(X);
      testCase.draw(X, y, idx);
    end
  end
end

