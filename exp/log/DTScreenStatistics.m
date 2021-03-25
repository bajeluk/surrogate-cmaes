classdef DTScreenStatistics < Observer
%SCREENSTATISTICS -- print statistics from DoubleTrainEC on screen without adaptivity
  properties
    verbosity
  end

  methods
    function obj = DTScreenStatistics(params)
      obj@Observer();
      verbosity = defopts(params, 'verbose', 5);
    end

    function notify(obj, ec, varargin)
    % get the interesting data from the generation, process and show them

      % get EC abbreviation
      switch class(ec)
        case 'DoubleTrainedEC'
          ecShort = 'DTS';
        case 'GenerationEC'
          ecShort = ' S ';
        case 'LinQuadEC'
          ecShort = 'lq';
        case 'LmmEC'
          ecShort = 'lmm';
        case 'ModelAssistedEC'
          ecShort = 'MA';
        otherwise
          ecShort = class(ec);
          if strcmp(ecShort(end-1:end), 'EC')
            ecShort = ecShort(1:end-2);
          end
      end
      % each 10 generations show the header
      if (mod(ec.cmaesState.countiter, 10) == 1)
      %           #####  iter /evals(or:p,b) | Dopt |rmseR | rnkR | rnk2 |rnkVal * | Mo nD nDiR |sigm^2| aErr |smooEr| orRat| aGain|
        fprintf('##%s## iter /evals(or:p,b) | D_fopt. | rmseRee | rnkR | rnk2 | .rankErrValid. | M  nData | sigma^2. | aErr |smooEr| orRat| aGain|\n', ...
                '#'*ones(1, numel(ecShort)));
      end
      model = '.';
      nTrainData = 0;
      % model successfully trained in current generation
      if (~isempty(ec.model) && ec.model.isTrained() ...
          && ec.model.trainGeneration == ec.cmaesState.countiter)
        model = '+';
        nTrainData = ec.model.getTrainsetSize();
      end
      % model successfully retrained in current generation
      if (~isempty(ec.retrainedModel) && ec.retrainedModel.isTrained() ...
          && ec.retrainedModel.trainGeneration == ec.cmaesState.countiter)
        model = '#';
        nTrainData = ec.retrainedModel.getTrainsetSize();
      end
      outputValues1 = [...
          ec.cmaesState.countiter, ... % generation number
          ec.counteval, ... % number of original-evaluated points so far
          sum(ec.pop.origEvaled), ... % number of original-evaluated points in current generation
          ec.nPresampledPoints, ... % number of presampled points
          ec.usedBestPoints, ...    % number of really original-evaluated points with the best predicted f-value
          ec.stats.fmin - ec.surrogateOpts.fopt, ... % difference between fitness minimum of original-evaluated points in current generation and original fitness minimum (?? verify ec.surrogateOpts.fopt)
          ec.stats.rmseReeval, ...     % RMSE of the re-evaluated point(s)
          ec.stats.rankErrReeval, ...  % rank error of population with re-evaluated point(s)
          ec.stats.rankErr2Models, ... % rank error between prediction of two models
          ec.stats.rankErrValid ];     % rank error between true fitness and model prediction on the validation set
      outputValues2 = [...
          nTrainData, ...                 % number of training points of the last successfully trained model (first training if second did not succeed)
          ec.stats.nDataInRange, ...      % number of all available data in the specified range from xmean using Mahalanobis distance (getDataNearPoint)
          ec.cmaesState.sigma^2, ...      % square of CMA-ES step-size
          ec.stats.adaptErr, ...          % last measured rankDiff during update()
          ec.stats.adaptSmoothedErr, ...  % smoothed error value used before fed into transfer function
          ec.stats.lastUsedOrigRatio, ... % OrigRatio param which was used (last) in the last generation
          ec.stats.adaptGain ];           % gain of original ratio (to be converted via min/max)
      outputValues1(isnan(outputValues1)) = 0.0;
      outputValues2(isnan(outputValues2)) = 0.0;
      %         #####  iter /evals(or,p) | Dopt |rmseR | rnkR | rnk2 |rnkVal * | Mo nD nDiR |sigm^2| aErr |smooEr| orRat| aGain|
      fprintf('=[%s]= %4d /%5d(%2d:%1d,%1d) | %.1e | %.1e | %.2f | %.2f | %.2f %s | %s %2d/%3d | %.2e | %.2f | %.2f | %.2f | %.2f |\n', ...
          ecShort, outputValues1(:), decorateKendall(1-2*ec.stats.rankErrValid), model, outputValues2(:) ...
          );

      % Log the utilization of the individ. models in the ModelPool
      if (isa(ec.model, 'ModelPool') && (mod(ec.cmaesState.countiter, 10) == 0))
        fprintf('Best Models History:  %s\n', num2str(ec.model.bestModelsHistory));
      end

    end
  end
end
