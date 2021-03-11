function tests = modelTestSetsTest
% test creating datasets from experiments using modelTestSets
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase) 
% test empty input

  % no input should not cause error
  modelTestSets()
  % empty input should return empty dataset
  verifyEqual(testCase, modelTestSets([]), {[]})
end

function testOutput(testCase) 
% test valid results

  % set experiment
  exp_id = 'exp_doubleEC_28_log_nonadapt';

  fun = 1;
  dim = 2;
  inst = 11;

  opts.datasetName = 'DTS_test';
  opts.isForData = true;
  opts.nSnapshotsPerRun = 25;
  opts.outputDirname = 'exp_deleteme';
  opts.rewriteResults = false;
  opts.sampleMethod = 'uniform_wor';

  ds = modelTestSets(exp_id, fun, dim, inst, opts);
  % check resulting dataset
  report = checkModelTestSets(ds);
  
  % verify emptyness of individual stats
  verifyEmpty(testCase, report.emptySet)
  verifyEmpty(testCase, report.testSetX.wrongSize)
  verifyEmpty(testCase, report.testSetX.duplicity)
  verifyEmpty(testCase, report.testSetY.wrongSize)
  
end