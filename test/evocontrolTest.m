function tests = evocontrolTest
  tests = functiontests(localfunctions);
end

function testSimpleCmaes(testCase)
  dim = 2;

  [xmin, fmin, counteval] = s_cmaes('fellii', 2*ones(dim,1), 2);

  verifyEqual(testCase, counteval, 300*dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-6);
end

function testNoModelNoEvoControl(testCase)
  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  dim = 2;

  [xmin, fmin, counteval] = s_cmaes('fellii', 2*ones(dim,1), 2, [], 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 300*dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-6);
end

function testGpModelGenEvoControl(testCase)
  dim = 2;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'rf';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = dim;     % 0..inf

  cmaesOpts = struct();
  % cmaesOpts.PopSize = 20;

  [xmin, fmin, counteval, stopflag] = s_cmaes('frosen', 2*ones(dim,1), 2, cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-7);
end
