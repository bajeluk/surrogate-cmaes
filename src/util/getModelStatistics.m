function surrogateStats = getModelStatistics(model, cmaesState, surrogateOpts, sampleOpts, counteval, varargin)
%GETMODELSTATISTICS prints and save the statistics about the currently trained model on testing data
% Input:
%       model   - model to be evaluated
%       cmaesState      - CMA-ES state variables (struct)
%       surrogateOpts   - options from experiment definition
%       sampleOpts      - sampling options
%       counteval       - number of evaluations so-far used
%       varargin{1}     - last ranking error in current population

  % CMA-ES state variables
  xmean = cmaesState.xmean;
  sigma = cmaesState.sigma;
  lambda = cmaesState.lambda;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  countiter = cmaesState.countiter;

  [~, xValidTest, ~] = ...
      sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
  surrogateStats = NaN(1,7);
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

    % Get actual errRankMuOnly()
    [~, sort1] = sort(yTest);
    ranking2   = ranking(yPredict);
    errRank = errRankMuOnly(ranking2(sort1), cmaesState.mu);

    [minstd minstdidx] = min(cmaesState.sigma*sqrt(cmaesState.diagC));
    [maxstd maxstdidx] = max(cmaesState.sigma*sqrt(cmaesState.diagC));
    surrogateStats = [rmse, kendall, cmaesState.sigma, ...
        max(cmaesState.diagD)/min(cmaesState.diagD), minstd, maxstd, errRank];
    varParams = cell2mat(varargin);
    surrogateStats((end+1):(end+length(varParams))) = varParams;
    
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
    % the numbers of evaluations which will trigger data saving:
    testingEvals = surrogateOpts.saveModelTrainingData;
    idxLastReached = find(counteval > testingEvals);
    if (~isempty(idxLastReached))
      idxLastReached = idxLastReached(end);
      evalsReached = surrogateOpts.saveModelTrainingData(idxLastReached);
      filename = sprintf([surrogateOpts.experimentPath filesep 'modeltrain_f%s_%d.mat'], surrogateOpts.expFileID, evalsReached);
      if (~exist(filename, 'file'))
        trainsetX = model.getDataset_X();
        trainsetY = model.getDataset_y();
        testsetX = xValidTest';
        testsetY = yTest;
        surrogateOpts.modelOpts.bbob_func = [];
        sampleOpts.xintobounds = [];
        save(filename, 'trainsetX', 'trainsetY', 'testsetX', 'testsetY', 'evalsReached', 'surrogateOpts', 'lambda', 'sigma', 'xmean', 'BD', 'diagD', 'kendall', 'rmse');
      end
    end
  end
end
