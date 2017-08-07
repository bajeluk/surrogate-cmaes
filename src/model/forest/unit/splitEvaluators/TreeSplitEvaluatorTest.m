classdef TreeSplitEvaluatorTest
  properties (Constant)
    drawEnabled = true
  end
  
  methods (Static)
    function draw(X, y, split)      
      if TreeSplitEvaluatorTest.drawEnabled
        stack = dbstack;
        figure('Name', stack(2).name);
        scatter3(X(split, 1), X(split, 2), y(split));
        hold on;
        scatter3(X(~split, 1), X(~split, 2), y(~split));
        hold off;
        legend('yLeftWanted', 'yRightWanted');
      end
    end
  end
end

