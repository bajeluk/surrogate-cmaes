classdef DoubleTrainedEC < EvolutionControl
  properties 
    model
    
    restrictedParam
    useDoubleTraining
  end
  
  methods 
    function obj = DoubleTrainedEC(surrogateOpts)
    % constructor
      obj.model = [];
      obj.restrictedParam = surrogateOpts.evoControlRestrictedParam;
      obj.useDoubleTraining = defopts(surrogateOpts, 'evoControlUseDoubleTraining', true);
    end
    
    function [fitness_raw, arx, arxvalid, arz, counteval, lambda, archive, surrogateStats, origEvaled] = runGeneration(obj, cmaesState, surrogateOpts, sampleOpts, archive, counteval, varargin)
    % Run one generation of double trained evolution control
      
      fitness_raw = [];
      arxvalid = [];
      surrogateStats = NaN(1, 10);
      origEvaled = false(1, cmaesState.lambda);

      % extract cmaes state variables
      xmean = cmaesState.xmean;
      sigma = cmaesState.sigma;
      lambda = cmaesState.lambda;
      BD = cmaesState.BD;
      dim = cmaesState.dim;
      mu = cmaesState.mu;
      countiter = cmaesState.countiter;
      
      % prepare the final population to be returned to CMA-ES
      finalPop = FinalPopulation(lambda, dim);
      
      obj.model = ModelFactory.createModel(surrogateOpts.modelType, surrogateOpts.modelOpts, xmean');

      if (isempty(obj.model))
        % model could not be created :( use the standard CMA-ES
        return;
      end
      
      minTrainSize = obj.model.getNTrainData();

      nArchivePoints = myeval(surrogateOpts.evoControlTrainNArchivePoints);
      [xTrain, yTrain] = archive.getDataNearPoint(nArchivePoints, ...
          xmean', surrogateOpts.evoControlTrainRange, sigma, BD);
      
      [fitness_raw, arx, arxvalid, arz, archive, counteval, xTrain, yTrain] = ...
        presample(minTrainSize, cmaesState, surrogateOpts, sampleOpts, archive, counteval, xTrain, yTrain, varargin{:});

      nPresampledPoints = size(arxvalid, 2);
      if (nPresampledPoints == lambda)
        origEvaled = true(1, lambda);
        return
      end
      finalPop = finalPop.addPoints(arxvalid, fitness_raw, arx, arz, nPresampledPoints);

      % train the model 
      obj.model = obj.model.train(xTrain, yTrain, cmaesState, sampleOpts);

      nLambdaRest = lambda - nPresampledPoints;
      
      if (~obj.model.isTrained())
        return
      end
      
      % sample new points
      [xExtend, xExtendValid, zExtend] = ...
          sampleCmaesNoFitness(sigma, nLambdaRest, cmaesState, sampleOpts);
      modelOutput = obj.model.getModelOutput(xExtendValid');
      % choose rho points with low confidence to reevaluate
      if any(strcmpi(obj.model.predictionType, {'sd2', 'poi', 'ei'}))
        % higher criterion is better (sd2, poi, ei)
        [~, pointID] = sort(modelOutput, 'descend');
      else
        % lower criterion is better (fvalues, lcb, fpoi, fei)
        [~, pointID] = sort(modelOutput, 'ascend');
      end
      reevalID = false(1, nLambdaRest);
      assert(obj.restrictedParam >= 0 && obj.restrictedParam <= 1, 'evoControlRestrictedParam out of bounds [0,1]');
      nReeval = ceil(nLambdaRest * obj.restrictedParam);
      reevalID(pointID(1:nReeval)) = true;
      xToReeval = xExtend(:, reevalID);
      xToReevalValid = xExtendValid(:, reevalID);
      zToReeval = zExtend(:, reevalID);

      % Debug -- for the correlation info into surrogateStats
      [yModel1, sd2Model1] = obj.model.predict(xExtendValid');
      nTrainPoints = length(yTrain);
      xReevalTrans = ( (sigma * BD) \ xToReevalValid);
      xDist = mean( sqrt(sum( xReevalTrans.^2 )) );

      % original-evaluate the chosen points
      [yNew, xNew, xNewValid, zNew, counteval] = ...
          sampleCmaesOnlyFitness(xToReeval, xToReevalValid, zToReeval, sigma, nReeval, counteval, cmaesState, sampleOpts, varargin{:});
      fprintf('counteval: %d\n', counteval)

      finalPop = finalPop.addPoints(xNewValid, yNew, xNew, zNew, nReeval);

      % update the Archive
      archive = archive.save(xNewValid', yNew', countiter);
      % the obj.models' dataset will be supplemented with this
      % new points during the next training using all the xTrain
      % calculate the models' precision
      yPredict = obj.model.predict(xNewValid');
      kendall = corr(yPredict, yNew', 'type', 'Kendall');
      rmse = sqrt(sum((yPredict' - yNew).^2))/length(yNew);
      fprintf('  model-gener.: %d preSamples, reevaluated %d pts, test RMSE = %f, Kendl. corr = %f.\n', nPresampledPoints, nReeval, rmse, kendall);
      surrogateStats(1:2) = [rmse, kendall];

      % TODO: restrictedParam adaptivity
%       alpha = surrogateOpts.evoControlAdaptivity;
%       if nReeval > 1
%         obj.restrictedParam = (1-alpha)*obj.restrictedParam + alpha*(1-kendall)/2;
%       else
%         obj.restrictedParam = (1-alpha)*obj.restrictedParam + alpha*rmse;
%       end
%       fprintf('Restricted param: %f\n', obj.restrictedParam);

      if ~all(reevalID)
        xTrain = [xTrain; xNewValid'];
        yTrain = [yTrain; yNew'];
        % train the model again
        retrainedModel = obj.model.train(xTrain, yTrain, cmaesState, sampleOpts);
        if (obj.useDoubleTraining && retrainedModel.isTrained())
          yNewRestricted = retrainedModel.predict((xExtendValid(:, ~reevalID))');

          % Debug -- for the correlation info into surrogateStats
          [yModel2, sd2Model2] = retrainedModel.predict(xExtendValid');
          [~, sort1] = sort(yModel1);
          ranking2   = ranking(yModel2);
          err = errRankMuOnly(ranking2(sort1), mu);
          fprintf('Rank error: %f\n', err);

          surrogateStats_ = getModelStatistics(retrainedModel, cmaesState, surrogateOpts, sampleOpts, counteval, err, nTrainPoints, xDist);
          surrogateStats(1:length(surrogateStats_)) = surrogateStats_;
        else
          % use values estimated by the old model
          fprintf('DoubleTrainedEC: The new model could (is not set to) be trained, using the not-retrained model.\n');
          yNewRestricted = yPredict(~reevalID);
          surrogateStats_ = getModelStatistics(obj.model, cmaesState, surrogateOpts, sampleOpts, counteval);
          surrogateStats(1:length(surrogateStats_)) = surrogateStats_;
        end

        finalPop = finalPop.addPoints(xExtendValid(:, ~reevalID), yNewRestricted, xExtend(:, ~reevalID), zExtend(:, ~reevalID), 0);

        yNew = [yNew, yNewRestricted'];
        xNew = [xNew, xExtend(:, ~reevalID)];
        xNewValid = [xNewValid, xExtendValid(:, ~reevalID)];
        zNew = [zNew, zExtend(:, ~reevalID)];
        
      end
              
      assert(finalPop.nPoints == lambda, 'There are not yet all lambda points prepared, but they should be!');

      % sort the returned solutions (the best to be first)
      [finalPop, popSortInd] = finalPop.sort;

      % shift the f-values:
      %   if the model predictions are better than the best original value
      %   in the model's dataset, shift ALL (!) function values
      %   Note: - all values have to be shifted in order to preserve predicted
      %           ordering of values
      %         - small constant is added because of the rounding errors
      %           when numbers of different orders of magnitude are summed
      fminDataset = min(obj.model.dataset.y);
      fminModel = finalPop.getMinModeled;
      diff = max(fminDataset - fminModel, 0);
      finalPop = finalPop.shiftY(1.000001*diff);

      % save the resulting re-evaluated population as the returning parameters
      fitness_raw = finalPop.y;
      arx = finalPop.arx;
      arxvalid = finalPop.x;
      arz = finalPop.arz;
      origEvaled = finalPop.origEvaled;
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
