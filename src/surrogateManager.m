function [fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, lambda] = surrogateManager(xmean, sigma, lambda, BD, diagD, countiter, fitfun_handle, inOpts, varargin)
% surrogateManager  controls sampling of new solutions and using a surrogate model
%
% @xmean, @sigma, @lambda, @BD, @diagD -- CMA-ES internal variables
% @countiter            the number of the current generation
% @fitfun_handle        handle to fitness function (CMA-ES uses string name of the function)
% @inOpts               options/settings for surrogate modelling
% @varargin             arguments for the fitness function (varargin from CMA-ES)
%
% returns:
% @surrogateStats       vector of numbers with variable model/surrogate statsitics

  persistent ec;
  % ec - one such instance for evolution control, persistent between
  % different calls of the surrogateManager() function

  persistent archive;           % archive of original-evaluated individuals

  % TODO: make an array with all the status variables from each generation

  % Defaults for surrogateOpts
  sDefaults.evoControl  = 'none';                 % none | individual | generation | doubletrained(restricted)
  sDefaults.sampleFcn   = @sampleCmaes;           % sampleCmaes | ??? TODO ???
  sDefaults.evoControlPreSampleSize       = 0.2;  % 0..1
  sDefaults.evoControlIndividualExtension = 20;   % 1..inf (reasonable 10-100)
  sDefaults.evoControlSamplePreprocessing = false;
  sDefaults.evoControlBestFromExtension   = 0.2;  % 0..1
  sDefaults.evoControlTrainRange          = 8;    % 1..inf (reasonable 1--20)
  sDefaults.evoControlSampleRange         = 1;    % 1..inf (reasonable 1--20)
  sDefaults.evoControlOrigGenerations     = 1;    % 1..inf
  sDefaults.evoControlModelGenerations    = 1;    % 0..inf
  sDefaults.evoControlTrainNArchivePoints = 0;
  sDefaults.evoControlValidatePoints      = 0;
  sDefaults.evoControlRestrictedParam     = 0.2;    % 0..1
  sDefaults.evoControlSwitchMode          = 'none'; % none | individual | generation | doubletrained(restricted)
  sDefaults.evoControlSwitchBound         = inf;    % 1 .. inf (reasonable 10--100)
  sDefaults.evoControlSwitchPopulation    = 1;      % 1 .. inf (reasonable 1--20)
  sDefaults.evoControlSwitchPopBound      = inf;    % 1 .. inf (reasonable 10--100)
  sDefaults.modelType = '';                         % gp | rf
  sDefaults.modelOpts = [];                         % model specific options

  surrogateStats = NaN(1, 2);

  % copy the defaults settings...
  surrogateOpts = sDefaults;
  % and replace those set in the surrogateOpts:
  for fname = fieldnames(inOpts)'
    surrogateOpts.(fname{1}) = inOpts.(fname{1});
  end

  assert(size(xmean,2) == 1, 'surrogateManager(): xmean is not a column vector!');
  dim = size(xmean,1);

  % switching evolution control
  if surrogateOpts.sampleOpts.counteval > surrogateOpts.evoControlSwitchBound*dim
    surrogateOpts.evoControl = surrogateOpts.evoControlSwitchMode;
    % EC type has changed -> create new instance of EvolutionControl
    ec = ECFactory.createEC(surrogateOpts);
  end
  
  % switching population size
  if surrogateOpts.sampleOpts.origPopSize == lambda && surrogateOpts.sampleOpts.counteval > surrogateOpts.evoControlSwitchPopBound*dim
    lambda = ceil(surrogateOpts.evoControlSwitchPopulation * lambda);
  end

  % evolution control -- use model? individual? generation? double-trained?
  if (strcmpi(surrogateOpts.evoControl, 'none'))
    % No model at all
    if (strcmpi(func2str(surrogateOpts.sampleFcn), 'samplecmaes'))
      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
    else
      error('surrogateManager: the only sampling method without model is "sampleCmaes()"');
    end
    return;
  end

  if (countiter == 1)
    archive = Archive(dim);
    ec = ECFactory.createEC(surrogateOpts);
  end
  
  cmaesState = struct( ...
    'xmean', xmean, ...
    'sigma', sigma, ...
    'lambda', lambda, ...
    'BD', BD, ...
    'diagD', diagD, ...
    'dim', dim, ...
    'fitfun_handle', fitfun_handle, ...
    'countiter', countiter, ...
    'sampleOpts', surrogateOpts.sampleOpts);
  
  [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = ec.runGeneration(cmaesState, surrogateOpts, archive, varargin);
  surrogateOpts.sampleOpts.counteval = counteval;
  
  if (size(fitness_raw, 2) < lambda)
    % the model was in fact not trained
    disp('surrogateManager(): the model was not successfully trained.');
    [yNew, xNew, xNewValid, zNew, counteval] = sampleCmaes(xmean, sigma, lambda - size(fitness_raw, 2), BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
    surrogateOpts.sampleOpts.counteval = counteval;
    archive = archive.save(xNewValid', yNew', countiter);

    % save the resulting re-evaluated population as the returning parameters
    fitness_raw = [fitness_raw yNew];
    arx = [arx xNew];
    arxvalid = [arxvalid xNewValid];
    arz = [arz zNew];
  end
  
  assert(min(fitness_raw) >= min(archive.y), 'Assertion failed: minimal predicted fitness < min in archive by %e', min(archive.y) - min(fitness_raw));

end