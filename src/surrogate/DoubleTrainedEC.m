classdef DoubleTrainedEC < EvolutionControl
  properties 
    model
    pop
    cmaesState
    
    restrictedParam
    useDoubleTraining
    retrainedModel
    stats
    archive
    counteval
    nPresampledPoints
    surrogateOpts
  end
  
  methods 
    function obj = DoubleTrainedEC(surrogateOpts, varargin)
    % constructor
      obj.model = [];
      obj.restrictedParam = defopts(surrogateOpts, 'evoControlRestrictedParam', 0.1);
      obj.useDoubleTraining = defopts(surrogateOpts, 'evoControlUseDoubleTraining', true);
      obj.pop = [];
      stats = struct( ...
          'fmin', Inf ...
          );
      obj.surrogateOpts = surrogateOpts;
    end

    function [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = finalizeGeneration(obj, sampleOpts, varargin)
      % fill the rest of the population with original evaluations
      [obj, fitness_raw, arx, arxvalid, arz, counteval] = ...
          obj.fillPopWithOrigFitness(sampleOpts, varargin);

      % calculate and print statistics
      % TODO: refactor all output into Observers
      surrogateStats = NaN(1,10);
      if (~isempty(obj.model) && obj.model.isTrained() ....
          && obj.model.trainGeneration == obj.cmaesState.countiter)
        % get ranking returned by the first model of points from which were
        % decided point for re-evaluation (i.e. phase != 0)
        afterPresample = (obj.pop.phase ~= 0);
        xAfterPresample = obj.pop.x(:,afterPresample);
        yModel1 = obj.model.predict(xAfterPresample');

        % calculate RMSE and possibly Kendall's coeff. of the re-evaluated point(s)
        % (phase == 1)
        xReeval = obj.pop.x(:,(obj.pop.phase == 1));
        yReeval = obj.pop.y(:,(obj.pop.phase == 1));
        yReevalModel = obj.model.predict(xReeval');
        rmse = sqrt(sum((yReevalModel' - yReeval).^2))/length(yReeval);
        kendall = corr(yReevalModel, yReeval', 'type', 'Kendall');
        fprintf('  model: %d preSamples, reevaluated %d pts, test RMSE = %f, Kendl. corr = %f.\n', obj.nPresampledPoints, length(yReeval), rmse, kendall);
        surrogateStats(1:2) = [rmse, kendall];

        % calculate distance from mean -- TODO: aren't they obj.pop.arz...?!
        xReevalTrans = ( (obj.cmaesState.sigma * obj.cmaesState.BD) \ xReeval);
        xDist = mean( sqrt(sum( xReevalTrans.^2 )) );

        % get ranking error between the first and the second model
        if (~isempty(obj.retrainedModel) && obj.retrainedModel.isTrained())
          yModel2 = obj.retrainedModel.predict(xAfterPresample');
          [~, sort1] = sort(yModel1);
          ranking2   = ranking(yModel2);
          obj.stats.rankErr = errRankMuOnly(ranking2(sort1), obj.cmaesState.mu);
          fprintf('Rank error: %f\n', obj.stats.rankErr);
          statModel = obj.retrainedModel;
        else
          % ranking error between two models cannot be calculated :(
          obj.stats.rankErr = NaN;
          statModel = obj.model;
        end

        % save (and output) different statistics
        surrogateStats_ = getModelStatistics(statModel, obj.cmaesState, obj.surrogateOpts, sampleOpts, obj.counteval, obj.stats.rankErr, statModel.getTrainsetSize(), xDist);
        surrogateStats(1:length(surrogateStats_)) = surrogateStats_;
      end
      origEvaled = obj.pop.origEvaled;
    end
    
    function [obj, fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control
      
      % extract cmaes state variables
      obj.cmaesState = cmaesState;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      dim = cmaesState.dim;
      BD = cmaesState.BD;
      obj.archive = archive;
      obj.counteval = counteval;
      obj.retrainedModel = [];
      
      % prepare the final population to be returned to CMA-ES
      obj.pop = Population(lambda, obj.cmaesState.dim);
      
      newModel = ModelFactory.createModel(obj.surrogateOpts.modelType, obj.surrogateOpts.modelOpts, obj.cmaesState.xmean');

      if (isempty(newModel))
        % model could not be created :(. Use the standard CMA-ES.
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = obj.finalizeGeneration(sampleOpts, varargin);
        return;
      end
      
      minTrainSize = newModel.getNTrainData();

      nArchivePoints = myeval(obj.surrogateOpts.evoControlTrainNArchivePoints);
      [xTrain, yTrain] = obj.archive.getDataNearPoint(nArchivePoints, ...
          obj.cmaesState.xmean', obj.surrogateOpts.evoControlTrainRange, sigma, BD);
      
      % Do pre-sample
      [y, arx, x, arz, ~, obj.counteval, xTrain, yTrain] = ...
        presample(minTrainSize, obj.cmaesState, obj.surrogateOpts, sampleOpts, obj.archive, obj.counteval, xTrain, yTrain, varargin{:});
      obj.nPresampledPoints = size(x, 2);
      % Add all points as original-evaluated
      obj.pop = obj.pop.addPoints(x, y, arx, arz, obj.nPresampledPoints, 0);        % Phase == 0

      if (obj.nPresampledPoints == lambda)
        % everything has been evaluated and saved, so return
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = obj.finalizeGeneration(sampleOpts, varargin);
        return;
      end

      % train the model 
      newModel = newModel.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
      if (~newModel.isTrained())
        % model cannot be trained :( -- return with orig-evaluated population
        [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = obj.finalizeGeneration(sampleOpts, varargin);
        % TODO: try the old model instead just orig-evaluating all the population
        return;
      else
        obj.model = newModel;
      end

      nLambdaRest = lambda - obj.nPresampledPoints;

      % sample new points
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(sigma, nLambdaRest, obj.cmaesState, sampleOpts);
      [modelOutput, yExtendModel] = obj.model.getModelOutput(xExtend');

      reevalID = obj.choosePointsForReevaluation(xExtend, modelOutput, yExtendModel);
      xToReeval = xExtendValid(:, reevalID);
      nToReeval = sum(reevalID);

      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, obj.counteval] = ...
          sampleCmaesOnlyFitness(xExtend(:, reevalID), xToReeval, zExtend(:, reevalID), sigma, nToReeval, obj.counteval, obj.cmaesState, sampleOpts, varargin{:});
      fprintf('counteval: %d\n', obj.counteval)

      obj.pop = obj.pop.addPoints(xNewValid, yNew, xNew, zNew, nToReeval, 1);   % Phase == 1

      % update the Archive
      obj.archive.save(xNewValid', yNew', obj.cmaesState.countiter);

      % TODO: restrictedParam adaptivity
%       alpha = obj.surrogateOpts.evoControlAdaptivity;
%       if nToReeval > 1
%         obj.restrictedParam = (1-alpha)*obj.restrictedParam + alpha*(1-kendall)/2;
%       else
%         obj.restrictedParam = (1-alpha)*obj.restrictedParam + alpha*rmse;
%       end
%       fprintf('Restricted param: %f\n', obj.restrictedParam);

      if ~all(reevalID)
        xTrain = [xTrain; xNewValid'];
        yTrain = [yTrain; yNew'];
        % train the model again
        obj.retrainedModel = obj.model.train(xTrain, yTrain, obj.cmaesState, sampleOpts);
        if (obj.useDoubleTraining && obj.retrainedModel.isTrained())
          yReeval = obj.retrainedModel.predict((xExtendValid(:, ~reevalID))');
        else
          fprintf('DoubleTrainedEC: The new model could (is not set to) be trained, using the not-retrained model.\n');
          yReeval = yExtendModel(~reevalID);
        end
        obj.pop = obj.pop.addPoints(xExtendValid(:, ~reevalID), yReeval, ...
            xExtend(:, ~reevalID), zExtend(:, ~reevalID), 0, 2);   % Phase == 2
      end
              
      assert(obj.pop.nPoints == lambda, 'There are not yet all lambda points prepared, but they should be!');

      % sort the returned solutions (the best to be first)
      [obj.pop, ~] = obj.pop.sort;

      % shift the f-values:
      %   if the model predictions are better than the best original value
      %   in the model's dataset, shift ALL (!) function values
      %   Note: - all values have to be shifted in order to preserve predicted
      %           ordering of values
      %         - small constant is added because of the rounding errors
      %           when numbers of different orders of magnitude are summed
      fminDataset = min(obj.model.dataset.y);
      fminModel = obj.pop.getMinModeled;
      diff = max(fminDataset - fminModel, 0);
      obj.pop = obj.pop.shiftY(1.000001*diff);

      % save the resulting re-evaluated population as the returning parameters
      [obj, fitness_raw, arx, arxvalid, arz, counteval, surrogateStats, origEvaled] = obj.finalizeGeneration(sampleOpts, varargin);
    end
    
    function [obj, yNew, xNew, xNewValid, zNew, counteval] = fillPopWithOrigFitness(obj, sampleOpts, varargin)
      %Fill the rest of the current population 'pop' (of class Population) with 
      % the original-evaluated individuals
      nToEval = obj.cmaesState.lambda - sum(obj.pop.nPoints);
      if (nToEval > 0)
        [y, arx, x, arz, obj.counteval] = sampleCmaes(obj.cmaesState, sampleOpts, ...
            nToEval, obj.counteval, varargin{:});
        obj.archive.save(x', y', obj.cmaesState.countiter);
        obj.pop = obj.pop.addPoints(x, y, arx, arz, nToEval, 3);        % Phase == 3
      end
      obj.pop.sort();
      yNew = obj.pop.y;
      xNewValid = obj.pop.x;
      xNew = obj.pop.arx;
      zNew = obj.pop.arz;
      counteval = obj.counteval;
    end

    function reevalID = choosePointsForReevaluation(obj, xExtend, modelOutput, yExtendModel)
    % choose points with low confidence to reevaluate
    %
    % returns:
    %   reevalID  --  bool vector which points from xExtend to reevaluate
      if any(strcmpi(obj.model.predictionType, {'sd2', 'poi', 'ei'}))
        % higher criterion is better (sd2, poi, ei)
        [~, pointID] = sort(modelOutput, 'descend');
      elseif (strcmpi(obj.model.predictionType, 'expectedrank'))
        % TODO it should work (yExtendModel, modelOutput) instead of re-predict the points again
        [y_m, sd2_m] = obj.model.predict(xExtend');
        pointID = expectedRankDiff(y_m, sd2_m, obj.cmaesState.mu, @errRankMuOnly);
        % Debug:
        % y_r = ranking(y_m);
        % fprintf('  Expected permutation of sorted f-values: %s\n', num2str(y_r(pointID)'));
      else
        % lower criterion is better (fvalues, lcb, fpoi, fei)
        [~, pointID] = sort(modelOutput, 'ascend');
      end
      nLambdaRest = size(xExtend, 2);
      reevalID = false(1, nLambdaRest);
      assert(obj.restrictedParam >= 0 && obj.restrictedParam <= 1, 'evoControlRestrictedParam out of bounds [0,1]');
      nToReeval = ceil(nLambdaRest * obj.restrictedParam);
      reevalID(pointID(1:nToReeval)) = true;
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
