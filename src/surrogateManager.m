function [fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, lambda, origEvaled] = surrogateManager(cmaesState, inOpts, sampleOpts, counteval, varargin)
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
% @lambda               pop. size (can be changed during this run
% @origEvaled           binary vector with true values on indices of individuals which
%                       are evaluated by the original fitness

  persistent ec;
  % ec - one such instance for evolution control, persistent between
  % different calls of the surrogateManager() function

  persistent archive;           % archive of original-evaluated individuals

  persistent observers;         % observers of EvolutionControls

  % TODO: make an array with all the status variables from each generation
  
  % CMA-ES state variables
  xmean = cmaesState.xmean;
  lambda = cmaesState.lambda;
  countiter = cmaesState.countiter;

  % Defaults for surrogateOpts
  sDefaults.evoControl  = 'none';                 % none | individual | generation | doubletrained(restricted)
  sDefaults.evoControlPreSampleSize       = 0.2;  % 0..1
  sDefaults.evoControlIndividualExtension = 20;   % 1..inf (reasonable 10-100)
  sDefaults.evoControlSamplePreprocessing = false;
  sDefaults.evoControlBestFromExtension   = 0.2;  % 0..1
  sDefaults.evoControlTrainRange          = 8;    % 1..inf (reasonable 1--20)
  sDefaults.evoControlSampleRange         = 1;    % 1..inf (reasonable 1--20)
  sDefaults.evoControlOrigGenerations     = 1;    % 1..inf
  sDefaults.evoControlModelGenerations    = 1;    % 0..inf
  sDefaults.evoControlTrainNArchivePoints = '15*dim';
  sDefaults.evoControlValidatePoints      = 0;
  sDefaults.evoControlRestrictedParam     = 0.2;    % 0..1
  sDefaults.evoControlAdaptivity          = 0.1;
  sDefaults.evoControlSwitchMode          = 'none'; % none | individual | generation | doubletrained(restricted)
  sDefaults.evoControlSwitchBound         = inf;    % 1 .. inf (reasonable 10--100)
  sDefaults.evoControlSwitchPopulation    = 1;      % 1 .. inf (reasonable 1--20)
  sDefaults.evoControlSwitchPopBound      = inf;    % 1 .. inf (reasonable 10--100)
  sDefaults.modelType = '';                         % gp | rf
  sDefaults.modelOpts = [];                         % model specific options

  surrogateStats = [];
  origEvaled = false(1, lambda);

  % copy the defaults settings...
  surrogateOpts = sDefaults;
  % and replace those set in the surrogateOpts:
  for fname = fieldnames(inOpts)'
    surrogateOpts.(fname{1}) = inOpts.(fname{1});
  end

  assert(size(xmean,2) == 1, 'surrogateManager(): xmean is not a column vector!');
  dim = size(xmean,1);
  cmaesState.dim = dim;

  % switching evolution control
  % TODO: consider removing this completely
  if counteval > surrogateOpts.evoControlSwitchBound*dim
    surrogateOpts.evoControl = surrogateOpts.evoControlSwitchMode;
    % EC type has changed -> create new instance of EvolutionControl
    % TODO: delete the old instances of 'ec' and 'observers'
    ec = ECFactory.createEC(surrogateOpts);
    [ec, observers] = ObserverFactory.createObservers(ec, surrogateOpts);
  end
  
  % switching population size
  if ((sampleOpts.origPopSize == lambda) && (counteval >= surrogateOpts.evoControlSwitchPopBound*dim))
    lambda = ceil(surrogateOpts.evoControlSwitchPopulation * lambda);
    cmaesState.lambda = lambda;
  end

  % construct Archive, EvolutionControl and its Observers
  % Note: independent restarts of the whole CMA-ES still clears the archive
  if (countiter == 1 && (isempty(ec) || (counteval < ec.counteval)))
    archive = Archive(dim);
    ec = ECFactory.createEC(surrogateOpts);
    [ec, observers] = ObserverFactory.createObservers(ec, surrogateOpts);
  end
  
  % run one generation according to evolution control
  [ec, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = ...
    ec.runGeneration(cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin{:});

  if (size(fitness_raw, 2) < lambda)
    % the model was in fact not trained
    % good EvolutionControl should not let come here!!!
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('EvolutionControl came back without full population of lambda points!');
    disp('It shouldn''t happen. Rest of points will be orig-evaluated.');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    [yNew, xNew, xNewValid, zNew, counteval] = sampleCmaes(cmaesState, sampleOpts, lambda - size(fitness_raw, 2), counteval, varargin{:});
    archive = archive.save(xNewValid', yNew', countiter);

    % save the resulting re-evaluated population as the returning parameters
    fitness_raw = [fitness_raw yNew];
    arx = [arx xNew];
    arxvalid = [arxvalid xNewValid];
    arz = [arz zNew];
    origEvaled((end-length(yNew)+1):end) = true;
  end

  assert(min(fitness_raw) >= min(archive.y), 'Assertion failed: minimal predicted fitness < min in archive by %e', min(archive.y) - min(fitness_raw));

  % check that the resulting points in arxvalid are inside bound
  % constraints
  inBounds = all(arxvalid == sampleOpts.xintobounds(...
      arxvalid, sampleOpts.lbounds, sampleOpts.ubounds));
  assert(all(inBounds), 'Assertion failed: arxvalid is out of bounds in generation %d', countiter);
end
