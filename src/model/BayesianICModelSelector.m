classdef BayesianICModelSelector < ICModelSelector
  %BAYESIANICMODELSELECTOR Extends ICModelSelector information criteria for
  %  Bayesian models.
  %
  %  Gelman et al., 2014, Bayesian Data Analysis, 3rd edition
    
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling

    % model selector-specific properties
    xMean
    models
    modelNames
    nModels
    bestIdx

    % ICModelSelector properties
    ic % information criterion to use
    modelsIC
  end

  methods (Access = protected)
    function calcICs(obj, generation)
      calcICs@ICModelSelector.calcICs(obj);

      obj.modelsIC.dic1(generation, :) = inf;
      obj.modelsIC.dic2(generation, :) = inf;
      obj.modelsIC.waic1(generation, :) = inf;
      obj.modelsIC.waic2(generation, :) = inf;
      obj.modelsIC.rhat(generation, :) = nan;

      for mdlIdx = 1:obj.nModels
        if ~obj.isTrained(generation, mdlIdx)
          continue;
        end

        mdl = obj.models{mdlIdx};

        % ---- DIC ----
        likEst = -mdl.getNegLogEst();
        likSample = mdl.getNegLogLPost();
        likSample = reshape(likSample, numel(likSample), 1);
        pd1 = 2 * (likEst - mean(likSample));
        pd2 = 2 * (var(likSample));

        dic1 = -2 * likEst + 2 * pd1;
        dic2 = -2 * likEst + 2 * pd2;
        
        obj.modelsIC.dic1(generation, mdlIdx) = dic1;
        obj.modelsIC.dic2(generation, mdlIdx) = dic2;

        % ---- WAIC ----
        [lppd, ppd] = mdl.getLogPredDens(obj.nSimuPost);

        % summations are over training cases, sample means and variances
        % are over posterior simulations
        pwaic1 = 2 * sum(lppd - mean(ppd, 2));
        pwaic2 = sum(var(ppd, 0, 2));
        lppd = sum(lppd);
        waic1 = -lppd + pwaic1;
        waic2 = -lppd + pwaic2;

        obj.modelsIC.waic1(generation, mdlIdx) = waic1;
        obj.modelsIC.waic2(generation, mdlIdx) = waic2;
        % TODO: Watanabe's Bayesian Information criterion
        % obj.modelsIC.wbic(generation, mdlIdx) = wbic;

        % ---- Gelman-Rubin convergence statistic ----
        chains = mdl.getChains();
        [n, m] = size(chains);
        if n > 0 && m > 1
          B = n * var(mean(chains, 1));
          W = mean(var(chains, 1));
          V = (n-1) * W / n + B / n;

          obj.modelsIC.rhat(generation, mdlIdx) = sqrt(V / W);
        end
      end
    end
  end
  
  methods (Access = public)
    function obj = BayesianICModelSelector(modelOptions, xMean)
      obj = obj@ICModelSelector(modelOptions, xMean);

      obj.ic = defopts(modelOptions, 'ic', 'waic2');
      obj.nSimuPost = defopts(modelOptions, 'nSimuPost', 100);
      obj.modelsIC.dic1 = zeros(1, obj.nModels);
      obj.modelsIC.dic2 = zeros(1, obj.nModels);
      obj.modelsIC.waic1 = zeros(1, obj.nModels);
      obj.modelsIC.waic2 = zeros(1, obj.nModels);
      obj.modelsIC.wbic = zeros(1, obj.nModels);
      obj.modelsIC.rhat = zeros(1, obj.nModels);
    end
  end
  
end

