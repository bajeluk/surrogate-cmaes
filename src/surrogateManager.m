function [fitness_raw, arx, arxvalid, arz, counteval] = surrogateManager(xmean, sigma, lambda, BD, diagD, countiter, fitfun_handle, fitargs, inOpts)
% surrogateManager  controls sampling of new solutions and using a surrogate model
%
% @xmean, @sigma, @lambda, @BD, @diagD -- CMA-ES internal variables
% @countiter            the number of the current generation
% @fitfun_handle        handle to fitness function (CMA-ES uses string name of the function)
% @fitargs              arguments for the fitness function (varargin from CMA-ES)
% @surrogateOpts        options/settings for surrogate modelling

  % Defaults for surrogateOpts
  sDefaults.evoControl  = 'none';               % none | individual | generation
  sDefaults.sampleFcn   = @sampleCmaes;         % sampleCmaes | ??? TODO ???

  % copy the defaults settings...
  surrogateOpts = sDefaults;
  % and replace those set in the surrogateOpts:
  for fname = fieldnames(inOpts)'
    surrogateOpts.(fname{1}) = inOpts.(fname{1});
  end

  % evolution control -- use model? individual? generation?
  if (strcmpi(surrogateOpts.evoControl, 'none'))
    % No model at all
    if (strcmpi(func2str(surrogateOpts.sampleFcn), 'samplecmaes'))
      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);
    else
      error('surrogateManager: the only sampling method without model is "sampleCmaes()"');
    end

  elseif (strcmpi(surrogateOpts.evoControl, 'individual'))
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

    % % surrogateOpts.evoControlOrigGenerations
    % % surrogateOpts.evoControlModelGenerations
    % cycleLength = surrogateOpts.evoControlOrigGenerations + surrogateOpts.evoControlModelGenerations;
    % 
    % if (generationEC.evaluateOriginal)
    %   [fitness_raw, arx, arxvalid, arz, counteval] = ...
    %       sampleCmaesOnlyFitness(arx, arxvalid, arz, xmean, sigma, lambda, BD, diagD, fitfun_handle, fitargs, surrogateOpts.sampleOpts);
    %   archive.save(arx, fitness_raw);
    %   gens = generationEC.getLastOriginalGenerationSession();
    %   [X, y] = archive.getDataFromGeneration(gens);

  else
    error('surrogateManager: wrong evolution control method');
  end
end
