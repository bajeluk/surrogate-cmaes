% MultiTrainedEC Test Class
%
% author: Lukas Bajer
% date: 2016-08-25
%
% MultiTrainedEC is a class with Multi-trained Evolution Control
% for the MTS-CMA-ES

classdef MultiTrainedECTest < matlab.unittest.TestCase

  methods(TestMethodSetup)
    % Functions which run before each test method
    
    function setupRandStream(testCase)
      % Setup reproducible results
      s = RandStream('mt19937ar', 'Seed', 3);
      RandStream.setGlobalStream(s);
    end
  end  
  
  methods (Test)
    % Indivudal test cases

    function testConstructor(testCase)
      surrogateOpts = struct();
      mtec = MultiTrainedEC(surrogateOpts);
    end
    
    function testMostProbableRankDiff_6_3(testCase)
      surrogateOpts = struct();
      mtec = MultiTrainedEC(surrogateOpts);

      n = 6; mu = 3;
      % yPredict = rand(1,n);
      % sd2Predict = 0.6*rand(1,n);
      yPredict =   [0.5508    0.70815  0.2909   0.51083  0.89295   0.89629]';
      sd2Predict = [0.075351  0.12435  0.03088  0.26449  0.017926  0.2741]';
      
      perm = mtec.mostProbableRankDiff(yPredict, sd2Predict, mu);
      % Expected error of [1] is 0.166667.
      % Expected error of [2] is 0.228913.
      % Expected error of [3] is 0.168245.
      % Expected error of [4] is 0.021866.
      % Expected error of [5] is 0.000000.
      % Expected error of [6] is 0.025076.
      yRnk = ranking(yPredict);
      testCase.verifyEqual(yRnk(perm), [2 3 1 6 4 5]');
    end
      
    function testMostProbableRankDiff_5_3(testCase)
      surrogateOpts = struct();
      mtec = MultiTrainedEC(surrogateOpts);

      n = 5; mu = 3;
      % yPredict = rand(1,n);
      % sd2Predict = 0.6*rand(1,n);
      yPredict =   [0.5508  0.7081  0.2909  0.5108  0.8929]';
      sd2Predict = [0.5378  0.0754  0.1243  0.0309  0.2645]';

      perm = mtec.mostProbableRankDiff(yPredict, sd2Predict, mu);
      % Expected error of [1] is 0.211396.
      % Expected error of [2] is 0.219554.
      % Expected error of [3] is 0.222688.
      % Expected error of [4] is 0.004119.
      % Expected error of [5] is 0.028147.
      yRnk = ranking(yPredict);
      testCase.verifyEqual(yRnk(perm), [3 2 1 5 4]');
    end

  end
end
