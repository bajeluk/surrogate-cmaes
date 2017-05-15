function [ok, fitness_raw, arx, arxvalid, arz, archive, counteval, xTrain, yTrain] = presample(minTrainSize, cmaesState, surrogateOpts, sampleOpts, archive, counteval, xTrain, yTrain, varargin)

  xmean = cmaesState.xmean;
  sigma = cmaesState.sigma;
  lambda = cmaesState.lambda;
  BD = cmaesState.BD;
  dim = cmaesState.dim;
  countiter = cmaesState.countiter;

  % The number of points to be 'pre-sampled'
  assert(surrogateOpts.evoControlPreSampleSize >= 0 && surrogateOpts.evoControlPreSampleSize <= 1, 'preSampleSize out of bounds [0,1]');
  maxPresampleSize = ceil(surrogateOpts.evoControlPreSampleSize * lambda);
  missingTrainSize = max(minTrainSize - size(xTrain, 1), 0);
  
  fitness_raw = []; arx = []; arxvalid = []; arz = [];

  if (missingTrainSize > maxPresampleSize)
    % TODO: shouldn't we use an old model?
    fprintf(2, '  Not enough data for training model.\n');
    ok = false;
    return;
  else
    ok = true;
  end

  if (missingTrainSize > 0)
    % pre-sample new points, preferably in areas where we don't have
    % the points yet
    expandedSigma = surrogateOpts.evoControlSampleRange * sigma;
    [~, arxvalid, arz] = ...
        sampleCmaesNoFitness(expandedSigma, dim*lambda, cmaesState, sampleOpts);
    [xPreSample, zPreSample] = SurrogateSelector.chooseDistantPoints(missingTrainSize, arxvalid', arz', xTrain, xmean, expandedSigma, BD);
    % evaluate the 'preSample' with the original fitness
    [fitness_raw, arx, arxvalid, arz, counteval] = ...
        sampleCmaesOnlyFitness(xPreSample, xPreSample, zPreSample, expandedSigma, missingTrainSize, counteval, cmaesState, sampleOpts, varargin{:});
    archive = archive.save(arxvalid', fitness_raw', countiter);
    xTrain = [xTrain; arxvalid'];
    yTrain = [yTrain; fitness_raw'];
    % the newModels' dataset will be supplemented with this
    % new points during the next training using all the xTrain
  end

end
