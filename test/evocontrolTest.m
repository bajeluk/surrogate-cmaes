function tests = evocontrolTest
  tests = functiontests(localfunctions);
end

function testNoModelEvoControl(testCase)
  [xmin, fmin, counteval] = s_cmaes('fellii', [2 2 2 2], 2);

  verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-10);
  verifyEqual(testCase, xmin, [0 0 0 0]', 'AbsTol', 1e-7);

  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  [xmin, fmin, counteval] = s_cmaes('fellii', [2 2 2 2], 2, [], 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-10);
  verifyEqual(testCase, xmin, [0 0 0 0]', 'AbsTol', 1e-7);

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.sampleFcn = @sampleCmaes;

  [xmin, fmin, counteval] = s_cmaes('fellii', [2 2 2 2], 2, [], 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-10);
  verifyEqual(testCase, xmin, [0 0 0 0]', 'AbsTol', 1e-7);
end