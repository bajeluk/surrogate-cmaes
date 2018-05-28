classdef BayesianICModelSelector < ICModelSelector
  %BAYESIANICMODELSELECTOR Extends ICModelSelector information criteria for
  %  Bayesian models.
  %
  %  Gelman et al., 2014, Bayesian Data Analysis, 3rd edition
  
  properties
    nSimuPost
  end

  methods (Access = protected)
    function obj = calcICs(obj, generation)
      obj = calcICs@ICModelSelector(obj, generation);

      obj.modelsIC.dic1(generation, :) = inf(1, obj.nModels);
      obj.modelsIC.dic2(generation, :) = inf(1, obj.nModels);
      obj.modelsIC.waic1(generation, :) = inf(1, obj.nModels);
      obj.modelsIC.waic2(generation, :) = inf(1, obj.nModels);
      obj.modelsIC.rhat(generation, :) = cell(1, obj.nModels);

      for mdlIdx = 1:obj.nModels
        if ~obj.modelIsTrained(generation, mdlIdx)
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
        if ~isempty(obj.nSimuPost)
          [lppd, ppd] = mdl.getLogPredDens(obj.nSimuPost);
        else
          [lppd, ppd] = mdl.getLogPredDens();
        end

        % summations are over training cases, sample means and variances
        % are over posterior simulations
        pwaic1 = 2 * sum(lppd - mean(ppd, 2));
        pwaic2 = sum(var(ppd, 0, 2));
        lppd = sum(lppd);
        waic1 = -lppd + pwaic1;
        waic2 = -lppd + pwaic2;

        obj.modelsIC.lppd(generation, mdlIdx) = -lppd;
        obj.modelsIC.waic1(generation, mdlIdx) = waic1;
        obj.modelsIC.waic2(generation, mdlIdx) = waic2;
        % TODO: Watanabe's Bayesian Information criterion
        obj.modelsIC.wbic(generation, mdlIdx) = 0;

        % ---- Gelman-Rubin convergence statistic ----
        chains = mdl.getChains();
        [ntot, d] = size(chains);
        M = mdl.mcmcNChains;
        assert(~rem(ntot, M));
        N = fix(ntot / M);

        if N > 0 && M > 1
          rhat = zeros(1, d);
          for j = 1:d
            % loop over estimands

            ch = zeros(N, M);
            for k = 1:M
              ch(:, k) = chains((k-1)*N+1:k*N, j);
            end

            % half each chain, throwing away remainders
            n = fix(N / 2);
            r = rem(N, 2);
            ch = [ch(1:n, :) ch(n+1:(N-r), :)];

            B = n * var(mean(ch, 1));
            W = mean(var(ch, 0, 1));
            V = (n-1) * W / n + B / n;
            rhat(j) = sqrt(V / W);
          end

          obj.modelsIC.rhat{generation, mdlIdx} = rhat;
        end
      end
    end
  end
  
  methods (Access = public)
    function obj = BayesianICModelSelector(modelOptions, xMean)
      obj = obj@ICModelSelector(modelOptions, xMean);

      obj.ic = defopts(modelOptions, 'ic', 'waic2');
      obj.nSimuPost = defopts(modelOptions, 'nSimuPost', []);
      obj.modelsIC.dic1 = zeros(1, obj.nModels);
      obj.modelsIC.dic2 = zeros(1, obj.nModels);
      obj.modelsIC.waic1 = zeros(1, obj.nModels);
      obj.modelsIC.waic2 = zeros(1, obj.nModels);
      obj.modelsIC.wbic = zeros(1, obj.nModels);
      obj.modelsIC.lppd = zeros(1, obj.nModels);
      obj.modelsIC.rhat = cell(1, obj.nModels);
    end
  end

end

