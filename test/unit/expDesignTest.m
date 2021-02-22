function tests = expDesignTest
% unit test for dOptimalExp function
  tests = functiontests(localfunctions);
end

function testEmptyInput(testCase)
  % test the expDesign with empty input

  % empty input should not cause error and return empty ids set
  verifyEmpty(testCase, expDesign());
end

function testDOptExistingInput(testCase)
  % test D-optimal design on existing experiment

  expName = 'exp_DTSmodels_meta_03_rf_validation';

  % test
  ids = expDesign(expName, 'dopt');
  verifyNotEmpty(testCase, ids)
end

function testLHSExistingInput(testCase)
  % test LHS design on existing experiment

  expName = 'exp_DTSmodels_meta_03_rf_validation';
  nCombs = 100;

  % test
  ids = expDesign(expName, 'lhs', 'nCombs', nCombs);
  [~, lastWarningId] = lastwarn;
  % return correct number of ids or raise warning
  verifyTrue(testCase, (numel(ids) == nCombs) || ...
    (strcmp(lastWarningId, 'scmaes:expDesign:notAchieveRes')))
end