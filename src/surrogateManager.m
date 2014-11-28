function [fitness_raw, arx, arxvalid, arz, counteval] = surrogateManager(xmean, sigma, lambda, BD, diagD, countiter, fitfun_handle, fitargs, inOpts)
% surrogateManager  controls sampling of new solutions and using a surrogate model
%
% @xmean, @sigma, @lambda, @BD, @diagD -- CMA-ES internal variables
% @countiter            the number of the current generation
% @fitfun_handle        handle to fitness function (CMA-ES uses string name of the function)
% @fitargs              arguments for the fitness function (varargin from CMA-ES)
% @surrogateOpts        options/settings for surrogate modelling

  persistent generationEC;              
  % generationEC - one such instance for evolution control, persistent between
  % different calls of the surrogateManager() function
  % TODO: test, if generationEC really works; GenerationEC is now derived from handle

  persistent lastModel;         % last successfully trained model

  % TODO: make an array with all the status variables from each generation

  % Defaults for surrogateOpts
  sDefaults.evoControl  = 'none';               % none | individual | generation
  sDefaults.sampleFcn   = @sampleCmaes;         % sampleCmaes | ??? TODO ???

  % copy the defaults settings...
  surrogateOpts = sDefaults;
  % and replace those set in the surrogateOpts:
  for fname = fieldnames(inOpts)'
    surrogateOpts.(fname{1}) = inOpts.(fname{1});
  end

  % Call this only once (during the first generation)
  if (countiter == 1)
    if (isfield(surrogateOpts, 'evoControlInitialGenerations'))
      initialGens = surrogateOpts.evoControlInitialGenerations;
    else
      initialGens = 0;
    end
    generationEC = GenerationEC(surrogateOpts.evoControlOrigGenerations, ...
      surrogateOpts.evoControlModelGenerations, initialGens);
    lastModel = [];
  end

  % evolution control -- use model? individual? generation?
  if (strcmpi(surrogateOpts.evoControl, 'none'))
    % No model at all
    if (strcmpi(func2str(surrogateOpts.sampleFcn), 'samplecmaes'))
      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);
    else
      error('surrogateManager: the only sampling method without model is "sampleCmaes()"');
    end
    return;
  end

  % some evolution control should be used (individual- or generation-based)
  newModel = ModelFactory.createModel(surrogateOpts.modelType, ...
      surrogateOpts.modelOpts, xmean);

  if (strcmpi(surrogateOpts.evoControl, 'individual'))
    % Individual-based evolution control
    % TODO: write this!
    warning('surrogateManager: Individual control is not yet written :(');
    warning('surrogateManager: Using "sampleCmaes()"...');
    [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);
    return;

  elseif (strcmpi(surrogateOpts.evoControl, 'generation'))
    % Generation-based evolution control
    [arx, arxvalid, arz] = ...
        sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);
    [fitness_raw, arx, arxvalid, arz, counteval] = ...
        sampleCmaesOnlyFitness(arx, arxvalid, arz, xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);

    if (generationEC.evaluateOriginal)
      [fitness_raw, arx, arxvalid, arz, counteval] = ...
          sampleCmaesOnlyFitness(arx, arxvalid, arz, xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);
      archive.save(arx, fitness_raw);
      if (~ generationEC.isNextOriginal())
        % we will switch to 'model'-mode in the next generation
        % prepare data for a new model
        nRequired = newModel.getNTrainData();
        X = []; y = [];
        if (length(fitness_raw) < nRequired)
          pointsWeNeed = nRequired - length(fitness_raw);
          gensWeNeed = ceil(pointsWeNeed / lambda);
          gens = generationEC.getLastOriginalGenerations(gensWeNeed);
          [X, y] = archive.getDataFromGeneration(gen);
        end
        X = [arx; X];
        y = [fitness_raw; y];
        if (length(y) >= nRequired)
          % we have got enough data for new model! hurraayh!
          newModel = newModel.train(X, y, xmean, countiter);
          % TODO: archive the lastModel...?
          lastModel = newModel;
        else
          % not enough training data :( -- continue with another
          % 'original'-evaluated generation
          generationEC = generationEC.holdOn();
        end
      end       % ~ generationEC.isNextOriginal()

    else
      % evalute the current population with the @lastModel
      % generationEC.evaluateModel() == true
      % TODO: implement other re-fitting strategies for an old model

      assert(~isempty(lastModel), 'surrogateManager(): we are asked to use an EMPTY MODEL!');
      
    end

  else
    error('surrogateManager(): wrong evolution control method');
  end
end
