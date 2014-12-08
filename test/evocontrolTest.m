function tests = evocontrolTest
  tests = functiontests(localfunctions);
end

function testSimpleCmaes(testCase)
  dim = 2;

  cmaesOpts.DispModulo = '5';
  cmaesOpts.StopFitness = 1e-8;
  [xmin, fmin, counteval] = s_cmaes('fellii', 2*ones(dim,1), 2, cmaesOpts);

  verifyEqual(testCase, counteval, 120*dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-3);
end

function testNoModelNoEvoControl(testCase)
  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  dim = 2;

  cmaesOpts.DispModulo = '5';
  cmaesOpts.StopFitness = 1e-8;
  [xmin, fmin, counteval] = s_cmaes('fellii', 2*ones(dim,1), 2, cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 120*dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-3);
end

function testGpModelGenEvoControl(testCase)
  dim = 2;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlOrigGenerations = 2;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = dim;     % 0..inf

  cmaesOpts = struct();
  % cmaesOpts.PopSize = 20;
  cmaesOpts.DispModulo = '5';
  cmaesOpts.StopFitness = 1e-8;

  [xmin, fmin, counteval, stopflag] = s_cmaes('frosen', 2*ones(dim,1), 2, cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(dim,1), 'AbsTol', 1e-7);
end
