classdef BaseTest < matlab.unittest.TestCase
  
  properties
    testFigure
  end
  
  
  methods (Access = protected)
    function [x] = fBase(testCase)
      x = 3;
    end
  end
  
  methods(TestMethodSetup)
    function createFigure(testCase)
        % comment
        testCase.testFigure = figure;
    end
  end

  methods(TestMethodTeardown)
    function closeFigure(testCase)
        close(testCase.testFigure);
    end
  end
  
  methods (Test)
    function testBase(testCase)
      rng('default');
      1+1
      testCase.fBase()
    end
  end
end