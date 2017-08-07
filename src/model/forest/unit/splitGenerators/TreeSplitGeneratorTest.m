classdef TreeSplitGeneratorTest < matlab.unittest.TestCase
  properties (Constant)
    drawEnabled = true
  end
  
  methods (Static)
    function [best, count, capturedVars] = findBest(X, y, generator)
      generator.reset(X, y);
      count = 0;
      best = struct('err', inf);
      while generator.hasNext()
        splitter = generator.next();
        idx = splitter(X);
        err = sum(y ~= idx);
        if err > length(y) / 2
          idx = ~idx;
          err = sum(y ~= idx);
        end
        if err < best.err
          best.err = err;
          best.splitter = splitter;
        end
        count = count + 1;
      end
      f = functions(best.splitter);
      capturedVars = f.workspace{1};
      
      if TreeSplitGeneratorTest.drawEnabled
        stack = dbstack;
        figure('Name', stack(2).name);
        scatter3(X(:, 1), X(:, 2), y);
        idx = best.splitter(X);
        hold on;
        XT = X;
        if isfield(capturedVars, 'W')
          XT = XT * capturedVars.W;
        end
        scatter3(XT(idx, 1), XT(idx, 2), y(idx));
        scatter3(XT(~idx, 1), XT(~idx, 2), y(~idx));
        hold off;
        legend('y', 'yLeft', 'yRight');
      end
    end
  end
end

