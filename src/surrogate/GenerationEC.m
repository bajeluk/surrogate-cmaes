classdef GenerationEC
  properties
    lastModel
    model
    
    origGenerations;
    modelGenerations;
    currentMode         = 'original';
    currentGeneration   = 1;
    lastOriginalGenerations = [];
    remaining           = 2;
  end

  methods
    function obj = GenerationEC(surrogateOpts)
      
      initialGens = defopts(surrogateOpts, 'evoControlInitialGenerations', 0);

      if (initialGens > 0)
        obj.currentMode = 'initial';
        obj.remaining = initialGens;
      else
        obj.currentMode = 'original';
        obj.remaining = surrogateOpts.evoControlOrigGenerations;
      end

      obj.origGenerations = surrogateOpts.evoControlOrigGenerations;
      obj.modelGenerations = surrogateOpts.evoControlModelGenerations;
      obj.currentGeneration   = 1;
      obj.lastModel = [];
      obj.model = [];
    end
    
    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats] = runGeneration(obj, cmaesState, surrogateOpts, archive, varargin)
      % Run one generation of individual evolution control

      fitness_raw = [];
      arx = [];
      arxvalid = [];
      arz = [];
      
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      BD = cmaesState.BD;
      diagD = cmaesState.diagD;
      fitfun_handle = cmaesState.fitfun_handle;
      countiter = cmaesState.countiter;
      counteval = surrogateOpts.sampleOpts.counteval;
      
      sampleSigma = surrogateOpts.evoControlSampleRange * sigma;

      if (obj.evaluateOriginal)
        %
        % original-evaluated generation
        %
        [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sampleSigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});

        surrogateOpts.sampleOpts.counteval = counteval;
        archive = archive.save(arxvalid', fitness_raw', countiter);
        if (~ obj.isNextOriginal())
          % we will switch to 'obj.model'-mode in the next generation
          % prepare data for a new model

          [obj.model, surrogateStats, isTrained] = trainGenerationECModel(obj.model, archive, xmean, sigma, lambda, BD, diagD, surrogateOpts, countiter);

          if (isTrained)
            % TODO: archive the obj.lastModel...?
            obj.lastModel = obj.model;
          else
            % not enough training data :( -- continue with another
            % 'original'-evaluated generation
            obj = obj.holdOn();
            return;
          end
        end       % ~ obj.isNextOriginal()

      else        % obj.evaluateModel() == true
        %
        % evalute the current population with the @obj.lastModel
        %
        % TODO: implement other re-fitting strategies for an old model

        if (isempty(obj.lastModel))
          warning('surrogateManager(): we are asked to use an EMPTY MODEL! Using CMA-ES.');
          [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sampleSigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
          surrogateOpts.sampleOpts.counteval = counteval;
          archive = archive.save(arxvalid', fitness_raw', countiter);
          return;
        end

        % generate the new population (to be evaluated by the model)
        [arx, arxvalid, arz] = ...
            sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);

        % generate validating population (for measuring error of the prediction)
        % this is with the *original* sigma
        [~, xValidValid, zValid] = ...
            sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);
        % shift the model (and possibly evaluate some new points newX, newY = f(newX) )
        % newX = []; newY = []; newZ = []; evals = 0;
        [shiftedModel, evals, newX, newY, newZ] = obj.lastModel.generationUpdate(xmean', xValidValid', zValid', surrogateOpts.evoControlValidatePoints, fitfun_handle, varargin{:});
        % count the original evaluations
        surrogateOpts.sampleOpts.counteval = surrogateOpts.sampleOpts.counteval + evals;
        counteval = surrogateOpts.sampleOpts.counteval;
        fitness_raw = zeros(1,lambda);
        % use the original-evaluated xValid points to the new generation:
        if (evals > 0)
          archive = archive.save(newX, newY, countiter);
          % calculate 'z' for the shifted archive near-mean point
          % because this near-mean is not sampled as Z ~ N(0,1)
          newZ(1,:) = ((BD \ (newX(1,:)' - xmean)) ./ sigma)';
          % save this point to the final population
          arx(:,1) = newX(1,:)';
          % this is a little hack :/ -- we suppose that all the newX are valid
          % but this should be true since 'newX' is derived from 'xValidValid'
          arxvalid(:,1) = newX(1,:)';
          arz(:,1) = newZ(1,:)';
          fitness_raw(1) = newY(1)';
          remainingIdx = 2:lambda;
        else
          remainingIdx = 1:lambda;
        end
        % calculate/predict the fitness of the not-so-far evaluated points
        if (~isempty(shiftedModel))
          % we've got a valid model, so we'll use it!
          [predict_fitness_raw, ~] = shiftedModel.predict(arx(:,remainingIdx)');
          fitness_raw(remainingIdx) = predict_fitness_raw';
          disp(['Model.generationUpdate(): We are using the model for ' num2str(length(remainingIdx)) ' individuals.']);
          % shift the predicted fitness: the best predicted fitness
          % could not be better than the so-far best fitness -- it would fool CMA-ES!
          % TODO: test if shifting ALL THE INDIVIDUALS (not only predicted) would help?
          bestFitnessArchive = min(archive.y);
          bestFitnessPopulation = min(fitness_raw);
          diff = max(bestFitnessArchive - bestFitnessPopulation, 0);
          fitness_raw = fitness_raw + diff;

          % DEBUG:
          fprintf('  test ');
          surrogateStats = getModelStatistics(shiftedModel, xmean, sigma, lambda, BD, diagD, surrogateOpts, countiter);

        else
          % we don't have a good model, so original fitness will be used
          [fitness_raw_, arx_, arxvalid_, arz_, counteval] = ...
              sampleCmaesOnlyFitness(arx(:,remainingIdx), arxvalid(:,remainingIdx), arz(:,remainingIdx), xmean, sigma, length(remainingIdx), BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
          surrogateOpts.sampleOpts.counteval = counteval;
          arx(:,remainingIdx) = arx_;
          arxvalid(:,remainingIdx) = arxvalid_;
          arz(:,remainingIdx) = arz_;
          fitness_raw(remainingIdx) = fitness_raw_;
          archive = archive.save(arxvalid_', fitness_raw_', countiter);

          % train a new model for the next generation
          [obj.model, surrogateStats, isTrained] = trainGenerationECModel(obj.model, archive, xmean, sigma, lambda, BD, diagD, surrogateOpts, countiter);

          if (isTrained)
            % TODO: archive the obj.lastModel...?
            obj.lastModel = obj.model;
            % leave the next generation as a model-evaluated:
            obj = obj.holdOn();
          else
            % not enough training data :( -- continue with
            % 'original'-evaluated generation
            obj = obj.setNextOriginal();
          end
        end
        % and set the next as original-evaluated (later .next() will be called)
      end
      obj = obj.next();

    end
    
    function [newModel, surrogateStats, isTrained] = trainGenerationECModel(model, archive, xmean, sigma, lambda, BD, diagD, surrogateOpts, countiter)
      
      surrogateStats = NaN(1, 2);
      % train the 'model' on the relevant data in 'archive'
      isTrained = false;

      trainSigma = surrogateOpts.evoControlTrainRange * sigma;
      nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);

      nRequired = model.getNTrainData();
      [X, y] = archive.getDataNearPoint(nArchivePoints, ...
          xmean', surrogateOpts.evoControlTrainRange, trainSigma, BD);
      if (length(y) >= nRequired)
        % we have got enough data for new model! hurraayh!
        sampleVariables = struct( ...
          'xmean', xmean, ...
          'sigma', sigma, ...
          'lambda', lambda, ...
          'BD', BD, ...
          'diagD', diagD, ...
          'sampleOpts', surrogateOpts.sampleOpts);
        % TODO: omit the unnecessary variables xmean, sigma and BD
        % as they are already in sampleVariables
        newModel = model.train(X, y, xmean', countiter, sigma, BD, sampleVariables);
        isTrained = (newModel.trainGeneration > 0);

        % DEBUG: print and save the statistics about the currently
        % trained model on testing data (RMSE and Kendall's correlation)
        if (isTrained)
          fprintf('  model trained on %d points, train ', length(y));
          surrogateStats = getModelStatistics(newModel, xmean, sigma, lambda, BD, diagD, surrogateOpts, countiter);
        end
      else
        newModel = model;
        isTrained = false;
      end
    end

    function result = evaluateOriginal(obj)
      % test whether evalute with the original function
      result = any(strcmp(obj.currentMode, {'original', 'initial'}));
    end

    function result = isNextOriginal(obj)
      % check whether there will be 'original' mode after calling next()
      result = (any(strcmp(obj.currentMode, {'original', 'initial'})) && (obj.remaining > 1)) ...
          || (strcmp(obj.currentMode, 'model') && obj.remaining == 1);
    end

    function result = evaluateModel(obj)
      % test whether evalute with a model
      result = strcmp(obj.currentMode, 'model');
    end 

    function obj = next(obj)
      % change the currentMode if all the generations from
      % the current mode have passed
      obj.remaining = obj.remaining - 1;
      switch obj.currentMode
        case 'initial'
          if (obj.remaining == 0)
            obj.currentMode = 'original';
            obj.remaining = obj.origGenerations;
          end
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
        case 'original'
          if (obj.remaining == 0)
            obj.currentMode = 'model';
            obj.remaining = obj.modelGenerations;
          end
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
        case 'model'
          if (obj.remaining == 0)
            obj.currentMode = 'original';
            obj.remaining = obj.origGenerations;
          end
        otherwise
          error('GenerationEC: wrong currentMode.');
      end
      obj.currentGeneration = obj.currentGeneration + 1;
    end

    function obj = holdOn(obj)
      % call this instead of next() if you want to
      % leave the current mode
      if (any(strcmp(obj.currentMode, {'original', 'initial'})))
        obj.lastOriginalGenerations = [obj.lastOriginalGenerations obj.currentGeneration];
      end
      obj.currentGeneration = obj.currentGeneration + 1;
    end

    function obj = setNextOriginal(obj)
      % set the next generation and currentMode to 'original'
      % later in the same generation, next() is expected to be called
      obj.currentMode = 'original';
      obj.remaining = 2;
    end

    function gens = getLastOriginalGenerations(obj, n)
      % get the numbers of the last n generations when the original
      % model was used
      startID = length(obj.lastOriginalGenerations) - n + 1;
      if (startID <= 0)
        disp('GenerationEC.getLastOriginalGenerations(): not enough data in the archive');
        startID = 1;
      end
      gens = obj.lastOriginalGenerations(startID:end);
    end
  end
end
