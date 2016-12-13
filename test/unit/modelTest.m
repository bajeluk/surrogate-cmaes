function tests = modelTest
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% initial settings

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  evoControlTest_fitness        = 'frosen';
  evoControlTest_dim            = 2;
  
  evoControlTest_cmaesOpts.DispModulo = '5';
  evoControlTest_cmaesOpts.StopFitness = 1e-8;
  evoControlTest_cmaesOpts.MaxFunEvals = 20*evoControlTest_dim;
  % original PopSize = '(4 + floor(3*log(N)))'
  evoControlTest_cmaesOpts.PopSize = (4 + floor(3*log(evoControlTest_dim)));
  
  cd(fullfile('..', '..'))
end

function testModelPredictionType(testCase)

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.modelOpts.predictionType = 'EI';
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = 2*evoControlTest_dim;     % 0..inf

  [xmin, fmin, counteval, stopflag] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(evoControlTest_dim,1), 'AbsTol', 1e-7);
end

function testModelTransformationType(testCase)

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'rf';
%   bbob_handlesF = benchmarks('handles');
%   surrogateOpts.modelOpts.bbob_func = bbob_handlesF{1};
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.modelOpts.predictionType = 'restricted';
  surrogateOpts.modelOpts.nTrees = 50;
  surrogateOpts.modelOpts.nBestPoints = 0;
  surrogateOpts.modelOpts.transformCoordinates = false;
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = 2*evoControlTest_dim;     % 0..inf

  [xmin, fmin, counteval, stopflag] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(evoControlTest_dim,1), 'AbsTol', 1e-7);
end