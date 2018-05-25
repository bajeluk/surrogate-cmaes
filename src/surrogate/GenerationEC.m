classdef GenerationEC < EvolutionControl & Observable
  properties
    model
    counteval
    
    lastModel
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
      obj.counteval = 0;
    end
    
    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
      % Run one generation of generation evolution control
      
      surrogateStats = NaN(1, 7);
      origEvaled = false(1, cmaesState.lambda);
      
      % extract cmaes state variables
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      BD = cmaesState.BD;
      fitfun_handle = cmaesState.fitfun_handle;
      countiter = cmaesState.countiter;
      obj.counteval = counteval;
      
      sampleSigma = surrogateOpts.evoControlSampleRange * sigma;
      
      obj.model = ModelFactory.createModel(surrogateOpts.modelType, surrogateOpts.modelOpts, xmean');

      if (obj.evaluateOriginal)
        %
        % original-evaluated generation
        %
        [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(cmaesState, sampleOpts, lambda, counteval, 'Archive', archive, varargin{:});
        origEvaled(1:end) = true;

        archive = archive.save(arxvalid', fitness_raw', countiter);
        if (~ obj.isNextOriginal())
          % we will switch to 'obj.model'-mode in the next generation
          % prepare data for a new model

          [obj, surrogateStats, isTrained] = obj.trainGenerationECModel(cmaesState, surrogateOpts, sampleOpts, archive, counteval);

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
          [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(cmaesState, sampleOpts, lambda, counteval, 'Archive', archive, varargin{:});
          archive = archive.save(arxvalid', fitness_raw', countiter);
          origEvaled(1:end) = true;
          return;
        end

        % generate the new population (to be evaluated by the model)
        [arx, arxvalid, arz] = ...
            sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
        origEvaled(1:end) = false;

        % generate validating population (for measuring error of the prediction)
        % this is with the *original* sigma
        [~, xValidValid, zValid] = ...
            sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
        % shift the model (and possibly evaluate some new points newX, newY = f(newX) )
        % newX = []; newY = []; newZ = []; evals = 0;
        [shiftedModel, evals, newX, newY, newZ] = obj.lastModel.generationUpdate(xmean', xValidValid', zValid', surrogateOpts.evoControlValidatePoints, fitfun_handle, varargin{:});
        % count the original evaluations
        counteval = counteval + evals;
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
          disp(['Model.generationUpdate(): We are using the model for ', num2str(length(remainingIdx)), ' individuals.']);
          
          % shift the f-values:
          %   if the model predictions are better than the best original value
          %   in the model's dataset, shift ALL (!) function values
          %   Note: - all values have to be shifted in order to preserve predicted
          %           ordering of values
          %         - small constant is added because of the rounding errors
          %           when numbers of different orders of magnitude are summed
          bestFitnessArchive = min(archive.y);
          bestFitnessPopulation = min(fitness_raw);
          diff = max(bestFitnessArchive - bestFitnessPopulation, 0);
          fitness_raw = fitness_raw + 1.000001*diff;

          % DEBUG:
          fprintf('  test ');
          surrogateStats = getModelStatistics(shiftedModel, cmaesState, surrogateOpts, sampleOpts, counteval);

        else
          % we don't have a good model, so original fitness will be used
          [fitness_raw_, arx_, arxvalid_, arz_, counteval] = ...
              sampleCmaesOnlyFitness(arx(:,remainingIdx), arxvalid(:,remainingIdx), arz(:,remainingIdx), sigma, length(remainingIdx), counteval, cmaesState, sampleOpts, 'Archive', archive, varargin{:});
          arx(:,remainingIdx) = arx_;
          arxvalid(:,remainingIdx) = arxvalid_;
          arz(:,remainingIdx) = arz_;
          fitness_raw(remainingIdx) = fitness_raw_;
          archive = archive.save(arxvalid_', fitness_raw_', countiter);

          % train a new model for the next generation
          [obj, surrogateStats, isTrained] = obj.trainGenerationECModel(cmaesState, surrogateOpts, sampleOpts, archive, counteval);

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
    
    function [obj, surrogateStats, isTrained] = trainGenerationECModel(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval)
      
      dim = cmaesState.dim;
      
      surrogateStats = NaN(1, 7);
      % train the 'model' on the relevant data in 'archive'
      isTrained = false;

      trainSigma = surrogateOpts.evoControlTrainRange * cmaesState.sigma;
      nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);

      nRequired = obj.model.getNTrainData();
      [X, y] = archive.getDataNearPoint(nArchivePoints, cmaesState.xmean', ...
        surrogateOpts.evoControlTrainRange, trainSigma, cmaesState.BD);
      if (length(y) >= nRequired)
        % we have got enough data for new model! hurraayh!
        obj.model = obj.model.train(X, y, cmaesState, sampleOpts);
        isTrained = (obj.model.trainGeneration > 0);

        % DEBUG: print and save the statistics about the currently
        % trained obj.model on testing data (RMSE and Kendall's correlation)
        if (isTrained)
          fprintf('  model trained on %d points, train ', length(y));
          surrogateStats = getModelStatistics(obj.model, cmaesState, surrogateOpts, sampleOpts, counteval);
        end
      else
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
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations, obj.currentGeneration];
        case 'original'
          if (obj.remaining == 0)
            obj.currentMode = 'model';
            obj.remaining = obj.modelGenerations;
          end
          obj.lastOriginalGenerations = [obj.lastOriginalGenerations, obj.currentGeneration];
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
        obj.lastOriginalGenerations = [obj.lastOriginalGenerations, obj.currentGeneration];
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

function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
