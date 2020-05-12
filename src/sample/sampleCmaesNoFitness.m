function [arx, arxvalid, arz] = sampleCmaesNoFitness(sigma, lambda, cmaesState, opts, modelMinimum)
% sampleCmaesNoFitness  generates lambda points according to CMA-ES scheme
%
% No point is evaluated by the fitness
% TODO: test this!

  % CMA-ES state variables
  xmean = cmaesState.xmean;
  BD = cmaesState.BD;
  diagD = cmaesState.diagD;
  
  % CMA-ES sampling options
  noiseReevals = opts.noiseReevals;
  bnd.isactive = opts.isBoundActive;
  lbounds = opts.lbounds;
  ubounds = opts.ubounds;
  flgDiagonalOnly = opts.flgDiagonalOnly;
  noiseHandling = opts.noiseHandling;
  xintobounds = opts.xintobounds;

  isBoundActive = any(lbounds > -Inf) || any(ubounds < Inf); 
  
  N = size(xmean, 1);
  % Generate and evaluate lambda offspring
  if exist('modelMinimum', 'var') && all(~isnan(modelMinimum))    
    modelMinimum = modelMinimum - cmaesState.xmean';
    rescaling_factor = 1;
    if mahalanobisNorm(modelMinimum, cmaesState.diagD, cmaesState.BD, cmaesState.sigma) > sqrt(N)
      rescaling_factor = normrnd(0, 1, [1, N]);
      rescaling_factor = rescaling_factor .^ 2;
      rescaling_factor = rescaling_factor(1);
      rescaling_factor = rescaling_factor .^ 0.5;
      rescaling_factor = rescaling_factor / mahalanobisNorm(modelMinimum, cmaesState.diagD, cmaesState.BD, cmaesState.sigma);
    end
        
    inject = modelMinimum * (rescaling_factor / cmaesState.sigma);
      
     
    arz = [randn(N, lambda-1) BD\inject'];
    
    %arz = randn(N, lambda);
    tmp = BD * arz;
    %tmp = [tmp(:,(1:end-1)), inject'];

    arx = repmat(xmean, 1, lambda) + sigma * tmp; % Eq. (1)   
  else
    arz = randn(N, lambda);
  

    if ~flgDiagonalOnly
        arx = repmat(xmean, 1, lambda) + sigma * (BD * arz); % Eq. (1)
    else
        arx = repmat(xmean, 1, lambda) + repmat(sigma * diagD, 1, lambda) .* arz; 
    end

    if noiseHandling 
        if noiseEpsilon == 0
            arx = [arx arx(:,1:noiseReevals)]; 
        elseif flgDiagonalOnly
            arx = [arx arx(:,1:noiseReevals) + ...
              repmat(noiseEpsilon * sigma * diagD, 1, noiseReevals) ...
              .* randn(N,noiseReevals)]; 
        else 
            arx = [arx arx(:,1:noiseReevals) + ...
              noiseEpsilon * sigma * ...
              (BD * randn(N,noiseReevals))]; 
        end
    end
  end
  % You may handle constraints here. You may either resample
  % arz(:,k) and/or multiply it with a factor between -1 and 1
  % (the latter will decrease the overall step size) and
  % recalculate arx accordingly. Do not change arx or arz in any
  % other way.

  if ~isBoundActive
    arxvalid = arx;
  else
    arxvalid = xintobounds(arx, lbounds, ubounds);
  end

end

function res = mahalanobisNorm(x, diagD, BD, sigma)
    D = diag(diagD);
    B = BD * inv(D);
    res = sum(((B' * x') ./ diagD) .^2 ) .^ 0.5;
    res = res / sigma;
end