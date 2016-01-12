function [fitness_raw, arx, arxvalid, arz, archive, counteval, xTrain, yTrain] = presample(minTrainSize, surrogateOpts, cmaesState, archive, xTrain, yTrain, varargin)

  xmean = cmaesState.xmean;
  sigma = cmaesState.sigma;
  lambda = cmaesState.lambda;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  dim = cmaesState.dim;
  fitfun_handle = cmaesState.fitfun_handle;
  countiter = cmaesState.countiter;
  counteval = surrogateOpts.sampleOpts.counteval;

  % The number of points to be 'pre-sampled'
  assert(surrogateOpts.evoControlPreSampleSize >= 0 && surrogateOpts.evoControlPreSampleSize <= 1, 'preSampleSize out of bounds [0,1]');
  maxPresampleSize = ceil(surrogateOpts.evoControlPreSampleSize * lambda);
  missingTrainSize = max(minTrainSize - size(xTrain, 1), 0);
  
  fitness_raw = []; arx = []; arxvalid = []; arz = [];

  if (missingTrainSize > maxPresampleSize)
    % TODO: shouldn't we use an old model?
    disp('surrogateManager(): not enough data for training model.');
    [fitness_raw, arx, arxvalid, arz, counteval] = sampleCmaes(xmean, sigma, lambda, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
    archive = archive.save(arxvalid', fitness_raw', countiter);
    return;
  end

  if (missingTrainSize > 0)
    % pre-sample new points, preferably in areas where we don't have
    % the points yet
    expandedSigma = surrogateOpts.evoControlSampleRange * sigma;
    [arx, ~, arz] = ...
        sampleCmaesNoFitness(xmean, expandedSigma, dim*lambda, BD, diagD, surrogateOpts.sampleOpts);
    [xPreSample, zPreSample] = SurrogateSelector.chooseDistantPoints(missingTrainSize, arx', arz', xTrain, xmean, expandedSigma, BD);
    % evaluate the 'preSample' with the original fitness
    [fitness_raw, arx, arxvalid, arz, counteval] = ...
        sampleCmaesOnlyFitness(xPreSample, xPreSample, zPreSample, xmean, expandedSigma, missingTrainSize, BD, diagD, fitfun_handle, surrogateOpts.sampleOpts, varargin{:});
    archive = archive.save(arxvalid', fitness_raw', countiter);
    xTrain = [xTrain; arxvalid'];
    yTrain = [yTrain; fitness_raw'];
    % the newModels' dataset will be supplemented with this
    % new points during the next training using all the xTrain
  end

end