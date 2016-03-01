function tests = alltest
  tests = functiontests(localfunctions);
end

function testStandardCmaesWorks(testCase)
  [xmin, fmin, counteval] = s_cmaes('fellii', [2 2 2 2], 2);

  verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-10);
  verifyEqual(testCase, xmin, [0 0 0 0]', 'AbsTol', 1e-7);

  [xmin, fmin, counteval] = s_cmaes('frosen', [2 2 2 2], 2, [], 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 2000, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-14);
  verifyEqual(testCase, xmin, [1 1 1 1]', 'AbsTol', 1e-4);
  
  [xmin, fmin, counteval] = s_cmaes('fellii', [2 2 2 2], 2, [], 'SurrogateOptions', surrogateOpts, 1000);
  
  verifyEqual(testCase, counteval, 1400, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-10);
  verifyEqual(testCase, xmin, [0 0 0 0]', 'AbsTol', 1e-7);
end
