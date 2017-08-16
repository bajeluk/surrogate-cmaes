classdef DerivedTest < BaseTest
  methods (Access = protected)
    function [x] = fDerived(testCase)
      x = 4;
    end
  end
  
  methods (Test)
    function testDerived(testCase)
      testCase.fBase()
      testCase.fDerived()
      testCase.testBase()
    end
  end
end