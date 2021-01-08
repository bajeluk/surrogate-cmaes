function tests = loadCompAlgTest
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase) 
% test empty input

  % no input should not cause error
  [algData, algNames, algColors] = loadCompAlg();
  % 6 algorithms should be on output
  verifySize(testCase, algData, [1, 6])
  verifySize(testCase, algNames, [1, 6])
  verifySize(testCase, algColors, [6, 3])
end

function testDefinedAlgInput(testCase)
% test defined algorithms input

  % input settings
  dataPath = fullfile('exp', 'pproc', 'compAlgMat.mat');
  funcSet.BBfunc = 1:24;
  funcSet.dims = [2, 3, 5, 10, 20];
  algs = {'cmaes', 'smac'};
  % run method
  [algData, algNames, algColors] = loadCompAlg(dataPath, funcSet, algs);
  
  % nAlgs algorithms should be on output
  nAlgs = numel(algs);
  verifySize(testCase, algData, [1, nAlgs])
  verifySize(testCase, algNames, [1, nAlgs])
  verifySize(testCase, algColors, [nAlgs, 3])
  
  % output should be as required
  verifyEqual(testCase, algNames, {'CMA-ES', 'SMAC'})
end