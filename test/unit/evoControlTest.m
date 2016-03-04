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
  evoControlTest_cmaesOpts.SaveVariables = false;
  evoControlTest_cmaesOpts.LogModulo = 0;
  
end

function teardownOnce(testCase)
% final commands
  
  % delete cmaes products
  if exist('variablescmaes.mat', 'file')
    delete('outcmaes*.dat')
    delete('variablescmaes.mat')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% no evolution control %
%%%%%%%%%%%%%%%%%%%%%%%%

function testSimpleCmaes(testCase)
% test cmaes part of s_cmaes

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testNoEvoControl(testCase)
% test evolution control switched off

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

%%%%%%%%%%%%%%%%%%%%%%%
% generation strategy %
%%%%%%%%%%%%%%%%%%%%%%%

function testGenEvoControl(testCase)
% test generation evolution control

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

%%%%%%%%%%%%%%%%%%%%%%%
% individual strategy %
%%%%%%%%%%%%%%%%%%%%%%%

function testIndEvoControl(testCase)
% test individual evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlIndividualExtension = 10;
  surrogateOpts.evoControlBestFromExtension = 0.1;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testIndEvoControlPreSample(testCase)
% test individual evolution control using presample

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0.2;
  surrogateOpts.evoControlIndividualExtension = 10;
  surrogateOpts.evoControlBestFromExtension = 0.1;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double-trained strategy %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testRestrEvoControl(testCase)
% test double-trained (restricted) evolution control with 'restricted' 
% settings - for compatibility reasons

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'restricted';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlRestrictedParam = 0.1;

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testDoubleTrainEvoControl(testCase)
% test double-trained evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'doubletrained';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.updaterParams = {'startRatio', 0.1};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testDoubleTrainEvoControlPresample(testCase)
% test double-trained evolution control using presample

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'doubletrained';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0.2;
  surrogateOpts.updaterParams = {'startRatio', 0.1};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testDoubleTrainECRMSEUpdate(testCase)
% test double-trained evolution control using rmse updater

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'doubletrained';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.updaterType = 'rmse';
  surrogateOpts.updaterParams = {'startRatio', 0.5};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

%%%%%%%%%%%%%%%%%%%%%%
% strategy switching %
%%%%%%%%%%%%%%%%%%%%%%

% TODO: recognize that switch was performed

function testSwitchGenToNoneEC(testCase)
% test switching from generation to none evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlOrigGenerations = 1;
  surrogateOpts.evoControlModelGenerations = 1;
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlSwitchMode = 'none';
  surrogateOpts.evoControlSwitchBound = round(0.5*evoControlTest_cmaesOpts.MaxFunEvals/evoControlTest_dim);

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testSwitchGenToDoubleTrainEC(testCase)
% test switching from generation to double-trained evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlOrigGenerations = 1;
  surrogateOpts.evoControlModelGenerations = 1;
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlSwitchMode = 'doubletrained';
  surrogateOpts.evoControlSwitchBound = round(0.5*evoControlTest_cmaesOpts.MaxFunEvals/evoControlTest_dim);
  surrogateOpts.updaterParams = {'startRatio', 0.1};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testSwitchDoubleTrainToNoneEC(testCase)
% test switching from double-trained to none evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'doubletrained';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlSwitchMode = 'none';
  surrogateOpts.evoControlSwitchBound = round(0.5*evoControlTest_cmaesOpts.MaxFunEvals/evoControlTest_dim);
  surrogateOpts.updaterParams = {'startRatio', 0.1};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population size switching %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testSwitchPopDoubleTrainEC(testCase)
% test switching from double-trained to none evolution control

  global evoControlTest_cmaesOpts;
  global evoControlTest_fitness;
  global evoControlTest_dim;

  surrogateOpts.evoControl = 'doubletrained';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.evoControlPreSampleSize = 0;
  surrogateOpts.evoControlSwitchPopulation = 2;
  surrogateOpts.evoControlSwitchPopBound = round(0.5*evoControlTest_cmaesOpts.MaxFunEvals/evoControlTest_dim);
  surrogateOpts.updaterParams = {'startRatio', 0.1};

  [xmin, fmin, counteval] = s_cmaes(evoControlTest_fitness, 2*ones(evoControlTest_dim,1), 2, evoControlTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [evoControlTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*evoControlTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end
%}
