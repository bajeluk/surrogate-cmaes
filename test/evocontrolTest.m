function tests = evocontrolTest
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  global evocontrolTest_cmaesOpts;
  global evocontrolTest_fitness;
  global evocontrolTest_dim;

  evocontrolTest_fitness        = 'frosen';
  evocontrolTest_dim            = 2;
  evocontrolTest_cmaesOpts.DispModulo = '5';
  evocontrolTest_cmaesOpts.StopFitness = 1e-8;

  % original PopSize = '(4 + floor(3*log(N)))'
  evocontrolTest_cmaesOpts.PopSize = (4 + floor(3*log(evocontrolTest_dim)));
end

function testSimpleCmaes(testCase)
  global evocontrolTest_cmaesOpts;
  global evocontrolTest_fitness;
  global evocontrolTest_dim;

  [xmin, fmin, counteval] = s_cmaes(evocontrolTest_fitness, 2*ones(evocontrolTest_dim,1), 2, evocontrolTest_cmaesOpts);

  verifyEqual(testCase, counteval, 300*evocontrolTest_dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  % verifyEqual(testCase, xmin, zeros(evocontrolTest_dim,1), 'AbsTol', 1e-3);
end

function testNoModelNoEvoControl(testCase)
  global evocontrolTest_cmaesOpts;
  global evocontrolTest_fitness;
  global evocontrolTest_dim;

  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  [xmin, fmin, counteval] = s_cmaes(evocontrolTest_fitness, 2*ones(evocontrolTest_dim,1), 2, evocontrolTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifyEqual(testCase, counteval, 300*evocontrolTest_dim, 'RelTol', 0.3);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-8);
  % verifyEqual(testCase, xmin, zeros(evocontrolTest_dim,1), 'AbsTol', 1e-3);
end

%{
function testGpModelGenEvoControl(testCase)
  global evocontrolTest_cmaesOpts;
  global evocontrolTest_fitness;
  global evocontrolTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = 2*evocontrolTest_dim;     % 0..inf

  [xmin, fmin, counteval, stopflag] = s_cmaes(evocontrolTest_fitness, 2*ones(evocontrolTest_dim,1), 2, evocontrolTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(evocontrolTest_dim,1), 'AbsTol', 1e-7);
end
%}

function testGpModelIndEvoControl(testCase)
  global evocontrolTest_cmaesOpts;
  global evocontrolTest_fitness;
  global evocontrolTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlPreSampleSize = 0.4;

  [xmin, fmin, counteval, stopflag] = s_cmaes(evocontrolTest_fitness, 2*ones(evocontrolTest_dim,1), 2, evocontrolTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(evocontrolTest_dim,1), 'AbsTol', 1e-7);
end
