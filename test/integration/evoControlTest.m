function tests = evoControlTest
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

function teardownOnce(testCase)
% final commands
  
  % delete cmaes products
  if exist('variablescmaes.mat', 'file')
    delete('outcmaes*.dat')
    delete('variablescmaes.mat')
  end
end

function testSimpleCmaes(testCase)
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testNoEvoControl(testCase)
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testGenEvoControl(testCase)
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = 2*evoControlTest_dim;     % 0..inf

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testIndEvoControl(testCase)
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlIndividualExtension = 10;
  surrogateOpts.evoControlBestFromExtension = 0.1;
  surrogateOpts.evoControlTrainRange = 4;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testIndEvoControlPreSample(testCase)
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0.2;
  surrogateOpts.evoControlIndividualExtension = 10;
  surrogateOpts.evoControlBestFromExtension = 0.1;
  surrogateOpts.evoControlTrainRange = 4;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testRestrEvoControl(testCase)
  % restricted - for compatibility reasons
  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'restricted';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlRestrictedParam = 0.1;
  surrogateOpts.evoControlTrainRange = 4;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end