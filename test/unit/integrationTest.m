function tests = integrationTest
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% initial settings

  global integrationTest_cmaesOpts;
  global integrationTest_fitness;
  global integrationTest_dim;

  integrationTest_fitness        = 'frosen';
  integrationTest_dim            = 2;
  
  integrationTest_cmaesOpts.DispModulo = '5';
  integrationTest_cmaesOpts.StopFitness = 1e-8;
  integrationTest_cmaesOpts.MaxFunEvals = 50*integrationTest_dim;
  % original PopSize = '(4 + floor(3*log(N)))'
  integrationTest_cmaesOpts.PopSize = (4 + floor(3*log(integrationTest_dim)));
  
  cd('..')
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
  global integrationTest_cmaesOpts;
  global integrationTest_fitness;
  global integrationTest_dim;

  [xmin, fmin, counteval] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts);

  verifySize(testCase, xmin, [integrationTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*integrationTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end

function testNoModelNoEvoControl(testCase)
  global integrationTest_cmaesOpts;
  global integrationTest_fitness;
  global integrationTest_dim;

  surrogateOpts.evoControl = 'none';
  surrogateOpts.sampleFcn = @sampleCmaes;

  [xmin, fmin, counteval] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  verifySize(testCase, xmin, [integrationTest_dim, 1])
  verifyLessThan(testCase, counteval, 1.5*integrationTest_cmaesOpts.MaxFunEvals)
  verifyNotEmpty(testCase, fmin);
end


%{
function testGpModelGenEvoControl(testCase)
  global integrationTest_cmaesOpts;
  global integrationTest_fitness;
  global integrationTest_dim;

  surrogateOpts.evoControl = 'generation';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
  surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
  surrogateOpts.evoControlValidatePoints = 2*integrationTest_dim;     % 0..inf

  [xmin, fmin, counteval, stopflag] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(integrationTest_dim,1), 'AbsTol', 1e-7);
end


 
function testGpModelIndEvoControl(testCase)
  global integrationTest_cmaesOpts;
  global integrationTest_fitness;
  global integrationTest_dim;

  surrogateOpts.evoControl = 'individual';
  surrogateOpts.modelType = 'gp';
  surrogateOpts.modelOpts.path = 'src/vendor/gpml-matlab-v3.2/';
  surrogateOpts.modelOpts.initScript = 'src/vendor/gpml-matlab-v3.2/startup.m';
  surrogateOpts.evoControlPreSampleSize = 0.4;
  surrogateOpts.evoControlIndividualExtension = 10;
  surrogateOpts.evoControlBestFromExtension = 0.1;
  surrogateOpts.evoControlTrainRange = 4;

  [xmin, fmin, counteval, stopflag] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);

  fprintf('\nCMA-ES ended:\n');
  celldisp(stopflag);

  % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
  verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
  % verifyEqual(testCase, xmin, zeros(integrationTest_dim,1), 'AbsTol', 1e-7);
end

%}

% function testModelPredictionType(testCase)
% 
%   global integrationTest_cmaesOpts;
%   global integrationTest_fitness;
%   global integrationTest_dim;
% 
%   surrogateOpts.evoControl = 'generation';
%   surrogateOpts.modelType = 'gp';
%   surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
%   surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
%   surrogateOpts.modelOpts.predictionType = 'EI';
%   surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
%   surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
%   surrogateOpts.evoControlValidatePoints = 2*integrationTest_dim;     % 0..inf
% 
%   [xmin, fmin, counteval, stopflag] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);
% 
%   fprintf('\nCMA-ES ended:\n');
%   celldisp(stopflag);
% 
%   % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
%   verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
%   % verifyEqual(testCase, xmin, zeros(integrationTest_dim,1), 'AbsTol', 1e-7);
% end

% function testModelTransformationType(testCase)
% 
%   global integrationTest_cmaesOpts;
%   global integrationTest_fitness;
%   global integrationTest_dim;
% 
%   surrogateOpts.evoControl = 'generation';
%   surrogateOpts.modelType = 'rf';
% %   bbob_handlesF = benchmarks('handles');
% %   surrogateOpts.modelOpts.bbob_func = bbob_handlesF{1};
%   surrogateOpts.modelOpts.path = '../gpeda/src/vendor/gpml-matlab-v3.2/';
%   surrogateOpts.modelOpts.initScript = '../gpeda/src/vendor/gpml-matlab-v3.2/startup.m';
%   surrogateOpts.modelOpts.predictionType = 'restricted';
%   surrogateOpts.modelOpts.nTrees = 50;
%   surrogateOpts.modelOpts.nBestPoints = 0;
%   surrogateOpts.modelOpts.transformCoordinates = false;
%   surrogateOpts.evoControlOrigGenerations = 1;      % 1..inf
%   surrogateOpts.evoControlModelGenerations = 1;     % 0..inf
%   surrogateOpts.evoControlValidatePoints = 2*integrationTest_dim;     % 0..inf
% 
%   [xmin, fmin, counteval, stopflag] = s_cmaes(integrationTest_fitness, 2*ones(integrationTest_dim,1), 2, integrationTest_cmaesOpts, 'SurrogateOptions', surrogateOpts);
% 
%   fprintf('\nCMA-ES ended:\n');
%   celldisp(stopflag);
% 
%   % verifyEqual(testCase, counteval, 1100, 'RelTol', 0.2);
%   verifyEqual(testCase, fmin, 0, 'AbsTol', 1e-7);
%   % verifyEqual(testCase, xmin, zeros(integrationTest_dim,1), 'AbsTol', 1e-7);
% end

