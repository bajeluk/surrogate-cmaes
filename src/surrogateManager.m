function [fitness_raw, arx, arxvalid, arz, counteval, surrogateStats] = surrogateManager(xmean, sigma, lambda, BD, diagD, countiter, fitfun_handle, inOpts, varargin)
% surrogateManager  controls sampling of new solutions and using a surrogate model
%
% @xmean, @sigma, @lambda, @BD, @diagD -- CMA-ES internal variables
% @countiter            the number of the current generation
% @fitfun_handle        handle to fitness function (CMA-ES uses string name of the function)
% @surrogateOpts        options/settings for surrogate modelling
% @varargin             arguments for the fitness function (varargin from CMA-ES)
%
% returns:
% @surrogateStats       vector of numbers with variable model/surrogate statsitics

  persistent generationEC;
  % generationEC - one such instance for evolution control, persistent between
  % different calls of the surrogateManager() function
  % TODO: test, if generationEC really works; GenerationEC is now derived from handle

  persistent lastModel;         % last successfully trained model
  persistent archive;           % archive of original-evaluated individuals

  % TODO: make an array with all the status variables from each generation

  % Defaults for surrogateOpts
  sDefaults.evoControl  = 'none';               % none | individual | generation
  sDefaults.sampleFcn   = @sampleCmaes;         % sampleCmaes | ??? TODO ???
  sDefaults.evoControlOrigGenerations = 1;      % 1..inf
  sDefaults.evoControlModelGenerations = 1;     % 0..inf
  sDefaults.modelType = '';                     % gp | rf
  sDefaults.modelOpts = [];                     % model specific options
  
  surrogateStats = NaN(1, 2);

  % copy the defaults settings...
  surrogateOpts = sDefaults;
  % and replace those set in the surrogateOpts:
  for fname = fieldnames(inOpts)'
    surrogateOpts.(fname{1}) = inOpts.(fname{1});
  end

  assert(size(xmean,2) == 1, 'surrogateManager(): xmean is not a column vector!');
  dim = size(xmean,1);

  % evolution control -- use model? individual? generation?
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

  % some evolution control should be used (individual- or generation-based)
  newModel = ModelFactory.createModel(surrogateOpts.modelType, ...
      surrogateOpts.modelOpts, xmean');
  if (countiter == 1)
    lastModel = [];
    archive = Archive(dim);
  end
  
  if (isempty(newModel))
    % model could not be created :( use the standard CMA-ES
    [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
    surrogateOpts.sampleOpts.counteval = counteval;
    archive = archive.save(arx', fitness_raw', countiter);
    return;
  end

  if (strcmpi(surrogateOpts.evoControl, 'individual'))
    % Individual-based evolution control
  
    nRequired = newModel.getNTrainData();
    nPreSample = ceil(surrogateOpts.evoControlPreSampleSize * lambda);
    fitness_raw = []; arx = []; arxvalid = []; arz = [];

    % TODO: let the model to be trained on larger dataset
    %       than is really needed if there are the data available
    if (nRequired > nPreSample)
      [xTrain, yTrain] = archive.getDataNearPoint((nRequired - nPreSample), ...
          xmean', surrogateOpts.evoControlTrainRange, sigma, BD);
    else
      xTrain = []; yTrain = [];
    end
    nToSample = nRequired - size(xTrain, 1);

    if (nToSample > nPreSample)
      % TODO: shouldn't we use an old model?
      disp('surrogateManager(): not enough data for training model.');
      [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
      archive = archive.save(arx', fitness_raw', countiter);
      return;
    end

    if (nToSample > 0)
      % pre-sample new points, preferably in areas where we don't have
      % the points yet
      [arx, ~, arz] = ...
          sampleCmaesNoFitness(xmean, sigma, 2*lambda, BD, diagD, surrogateOpts.sampleOpts);
      [xPreSample, zPreSample] = SurrogateSelector.chooseDistantPoints(nToSample, arx', arz', xTrain, xmean, sigma, BD);
      % evaluate the 'preSample' with the original fitness
      [fitness_raw, arx, arxvalid, arz, counteval] = ...
          sampleCmaesOnlyFitness(xPreSample, xPreSample, zPreSample, xmean, sigma, nToSample, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
      archive = archive.save(arx', fitness_raw', countiter);
      xTrain = [xTrain; arx'];
      yTrain = [yTrain; fitness_raw'];
    end
    % train the model
    newModel = newModel.train(xTrain, yTrain, xmean', countiter);
    % sample the enlarged population of size 'gamma * (lambda - nToSample)'
    extendSize = ceil(surrogateOpts.evoControlIndividualExtension ...
        * (lambda - nToSample));
    [xExtend, xExtendValid, zExtend] = ...
        sampleCmaesNoFitness(xmean, sigma, extendSize, BD, diagD, surrogateOpts.sampleOpts);
    % calculate the model prediction for the extended population
    yExtend = newModel.predict(xExtend');

    nBest = min(ceil(lambda*surrogateOpts.evoControlBestFromExtension), lambda - nToSample - 1);
    nCluster = lambda - nToSample - nBest;
    [xToReeval, xToReevalValid, zToReeval] = ...
        SurrogateSelector.choosePointsToReevaluate(...
        xExtend, xExtendValid, zExtend, yExtend, nBest, nCluster);

    % original-evaluate the chosen points
    [yNew, xNew, xNewValid, zNew, counteval] = ...
        sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, xmean, sigma, nCluster + nBest, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
    surrogateOpts.sampleOpts.counteval = counteval;
    archive = archive.save(xNew', yNew', countiter);
    yPredict = newModel.predict(xNew');
    kendall = corr(yPredict, yNew', 'type', 'Kendall');
    rmse = sqrt(sum((yPredict' - yNew).^2))/length(yNew);
    fprintf('  model: %d preSamples, reevaluated %d pts, RMSE = %f, Kendl. corr = %f.\n', nToSample, nCluster + nBest, rmse, kendall);
    surrogateStats = [rmse kendall];
    % TODO: control the evolution process according to the model precision

    % save the resulting re-evaluated population as the returning parameters
    fitness_raw = [fitness_raw yNew];
    arx = [arx xNew];
    arxvalid = [arxvalid xNewValid];
    arz = [arz zNew];

  elseif (strcmpi(surrogateOpts.evoControl, 'generation'))
    % Generation-based evolution control

    % Call this only once (during the first generation)
    if (countiter == 1)
      if (isfield(surrogateOpts, 'evoControlInitialGenerations'))
        initialGens = surrogateOpts.evoControlInitialGenerations;
      else
        initialGens = 0;
      end
      generationEC = GenerationEC(surrogateOpts.evoControlOrigGenerations, ...
        surrogateOpts.evoControlModelGenerations, initialGens);
    end

    [arx, arxvalid, arz] = ...
        sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);

    if (generationEC.evaluateOriginal)
      [fitness_raw, arx, arxvalid, arz, counteval] = ...
          sampleCmaesOnlyFitness(arx, arxvalid, arz, xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
      surrogateOpts.sampleOpts.counteval = counteval;
      archive = archive.save(arx', fitness_raw', countiter);
      if (~ generationEC.isNextOriginal())
        % we will switch to 'model'-mode in the next generation
        % prepare data for a new model
        nRequired = newModel.getNTrainData();
        X = []; y = [];
        if (length(fitness_raw) < nRequired)
          [X, y] = archive.getDataNearPoint(nRequired - length(fitness_raw), ...
              xmean', surrogateOpts.evoControlTrainRange, sigma, BD);
        end
        X = [arx'; X];
        y = [fitness_raw'; y];
        if (length(y) >= nRequired)
          % we have got enough data for new model! hurraayh!
          newModel = newModel.train(X, y, xmean', countiter);
          % TODO: archive the lastModel...?
          lastModel = newModel;
        else
          % not enough training data :( -- continue with another
          % 'original'-evaluated generation
          generationEC = generationEC.holdOn();
          return;
        end
      end       % ~ generationEC.isNextOriginal()

    else        % generationEC.evaluateModel() == true
      % evalute the current population with the @lastModel
      % 
      % TODO: implement other re-fitting strategies for an old model

      if (isempty(lastModel))
        warning('surrogateManager(): we are asked to use an EMPTY MODEL! Using CMA-ES.');
        [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
        surrogateOpts.sampleOpts.counteval = counteval;
        archive = archive.save(arx', fitness_raw', countiter);
        return;
      end

      % generate validating population (for measuring error of the prediction)
      [xValid, ~, zValid] = ...
          sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);
      % shift the model (and possibly evaluate some new points newX, newY = f(newX) )
      % newX = []; newY = []; newZ = []; evals = 0;
      [shiftedModel, evals, newX, newY, newZ] = lastModel.shiftReevaluate(xmean', xValid', zValid', surrogateOpts.evoControlValidatePoints, fitfun_handle, varargin{:});
      % count the original evaluations
      surrogateOpts.sampleOpts.counteval = surrogateOpts.sampleOpts.counteval + evals;
      counteval = surrogateOpts.sampleOpts.counteval;
      % use the original-evaluated xValid points to the new generation:
      xValidUsedIdx = 1:(evals-1);
      if (evals > 0)
        archive = archive.save(newX, newY, countiter);
        % calculate 'z' for the shifted archive near-mean point
        newZ(1,:) = ((BD \ (newX(1,:)' - xmean)) ./ sigma)';
        if (~isempty(xValidUsedIdx))
          % some of the 'xValid' points were evaluated, too, so save them
          % to the final population 'arx'
          arx(:,xValidUsedIdx) = newX(2:end,:)';
          % this is a hack :/ -- we suppose that all the newX are valid
          arxvalid(:,xValidUsedIdx) = newX(2:end,:)';
          arz(:,xValidUsedIdx) = newZ(2:end,:)';
          fitness_raw(xValidUsedIdx) = newY(2:end)';
        end
      end
      % calculate/predict the fitness of the not-so-far evaluated points
      remainingIdx = 1:lambda;
      remainingIdx(xValidUsedIdx) = [];
      if (~isempty(shiftedModel))
        % we've got a valid model, so we'll use it!
        [fitness_raw_, ~] = shiftedModel.predict(arx(:,remainingIdx)');
        disp(['Model.shiftReevaluate(): We are using the model for ' num2str(length(remainingIdx)) ' individuals.']);
        % shift the predicted fitness: the best predicted fitness
        % could not be better than the so-far best fitness -- it would fool CMA-ES!
        % TODO: test if shifting ALL THE INDIVIDUALS (not only predicted) would help?
        bestFitnessDataset = min(archive.y);
        bestFitnessPredicted = min(fitness_raw_);
        diff = max(bestFitnessDataset - bestFitnessPredicted, 0);
        fitness_raw_ = fitness_raw_' + diff;
      else
        % we don't have a model, so original fitness will be used
        [fitness_raw_, arx_, arxvalid_, arz_, counteval] = ...
            sampleCmaesOnlyFitness(arx(:,remainingIdx), arxvalid(:,remainingIdx), arz(:,remainingIdx), xmean, sigma, length(remainingIdx), BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
        arx(:,remainingIdx) = arx_;
        arxvalid(:,remainingIdx) = arxvalid_;
        arz(:,remainingIdx) = arz_;
        generationEC = generationEC.setNextOriginal();
      end
      fitness_raw(remainingIdx) = fitness_raw_;
      % and set the next as original-evaluated (later .next() will be called)
    end
    generationEC = generationEC.next();
    %
    % end of Generation-based evolution control

    % DEBUG
    % disp(fitness_raw);
    
  else
    error('surrogateManager(): wrong evolution control method');
  end
end
