function surrogateStats = getModelStatistics(model, cmaesState, surrogateOpts)
% print and save the statistics about the currently
% trained model on testing data

  xmean = cmaesState.xmean;
  sigma = cmaesState.sigma;
  lambda = cmaesState.lambda;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  countiter = cmaesState.countiter;

  [~, xValidTest, ~] = ...
      sampleCmaesNoFitness(xmean, sigma, lambda, BD, diagD, surrogateOpts.sampleOpts);
  surrogateStats = [NaN NaN];
  if (isfield(surrogateOpts.modelOpts, 'bbob_func'))
    preciseModel = ModelFactory.createModel('bbob', surrogateOpts.modelOpts, xmean');
    yTest = preciseModel.predict(xValidTest');
    yPredict = model.predict(xValidTest');
    kendall = corr(yPredict, yTest, 'type', 'Kendall');
    rmse = sqrt(sum((yPredict - yTest).^2))/length(yPredict);
    fprintf('  test RMSE = %f, Kendl. corr = %f. ', rmse, kendall);

    % decorate the kendall rho coefficient :)
    kendallInStars = floor(abs((kendall) * 5));
    if (kendallInStars == 0)
      stars = '[   o   ]';
    elseif (isnan(kendallInStars))
      stars = '[! NaN !]';
    else
      if (kendall > 0) 
        mark = '*'; 
      else
        mark = '-';
      end
      space = ' ';
      stars = sprintf('[ %s%s ]', mark(ones(1,kendallInStars)), space(ones(1,5-kendallInStars)));
    end
    fprintf('%s\n', stars);

    surrogateStats = [rmse kendall];
    
    % experimental
    % coef = sqrt(sum(((yPredict - yTest)/norm(yTest)).^2))/length(yPredict);
    % fprintf('\ncoef: %f\n\n', coef);
  else
    fprintf('\n');
  end

  % save the training and testing data for model-training enhancements
  % if ... the model is fresh
  %    ... and we'd like to save the training data
  if (model.trainGeneration == (countiter - 1) ...
      && isfield(surrogateOpts, 'saveModelTrainingData') ...
      && isfield(surrogateOpts, 'experimentPath') ...
      && surrogateOpts.saveModelTrainingData)
    currentEvals = surrogateOpts.sampleOpts.counteval;
    % the numbers of evaluations which will trigger data saving:
    testingEvals = surrogateOpts.saveModelTrainingData;
    idxLastReached = find(currentEvals > testingEvals);
    if (~isempty(idxLastReached))
      idxLastReached = idxLastReached(end);
      evalsReached = surrogateOpts.saveModelTrainingData(idxLastReached);
      filename = sprintf([surrogateOpts.experimentPath filesep 'modeltrain_f%s_%d.mat'], surrogateOpts.expFileID, evalsReached);
      if (~exist(filename, 'file'))
        trainsetX = model.dataset.X;
        trainsetY = model.dataset.y;
        testsetX = xValidTest';
        testsetY = yTest;
        surrogateOpts.modelOpts.bbob_func = [];
        surrogateOpts.sampleOpts.xintobounds = [];
        save(filename, 'trainsetX', 'trainsetY', 'testsetX', 'testsetY', 'evalsReached', 'surrogateOpts', 'lambda', 'sigma', 'xmean', 'BD', 'diagD', 'kendall', 'rmse');
      end
    end
  end
end