function [newstate, break_result] = cmaes_finalize(state)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE
[fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval,countevalNaN,irun ...
    flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals,stopMaxIter,stopFunEvals,stopIter, ... 
    stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues, ... 
    arrEqualFunvals, flgDiagonalOnly, flgActiveCMA,noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo, ... 
    flgscience,flgsaving,strsaving,flgsavingfinal, savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B,BD, ... 
    C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1,ccovmu,ccov1_sep, ccovmu_sep, ... 
    damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries,noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg, ... 
    stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop] = cmaes_loadState(state);
break_result = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE


% Evaluate xmean and return best recent point in xmin
fmin = fitness.raw(1);
xmin = arxvalid(:, fitness.idx(1)); % Return best point of last generation.
if length(stopflag) > sum(strcmp(stopflag, 'stoptoresume')) % final stopping
  out.solutions.mean.f = ...
      feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
  counteval = counteval + 1;
  out.solutions.mean.evals = counteval;
  if out.solutions.mean.f < fitness.raw(1)
    fmin = out.solutions.mean.f;
    xmin = xintobounds(xmean, lbounds, ubounds); % Return xmean as best point
  end
  if out.solutions.mean.f < out.solutions.bestever.f
    out.solutions.bestever = out.solutions.mean; % Return xmean as bestever point
    out.solutions.bestever.x = xintobounds(xmean, lbounds, ubounds); 
    bestever = out.solutions.bestever;
  end
end

% Save everything and display final message
if flgsavingfinal
  clear idx; % prevents error under octave
  if ~isempty(strsaving) && ~isoctave
    save('-mat', strsaving, opts.SaveFilename); % for inspection and possible restart	
  else 
    save('-mat', opts.SaveFilename);    % for inspection and possible restart
  end
  message = [' (saved to ' opts.SaveFilename ')'];
else
  message = [];
end

if flgdisplay
  disp(['#Fevals:   f(returned x)   |    bestever.f     | stopflag' ...
        message]);
  if isoctave
    strstop = stopflag(:); 
  else
      strcat(stopflag(:), '.');
  end
  strstop = stopflag(:); %strcat(stopflag(:), '.');
  disp([repmat(' ',1,6-floor(log10(counteval))) ...
        num2str(counteval, '%6.0f') ': ' num2str(fmin, '%.11e') ' | ' ...
        num2str(out.solutions.bestever.f, '%.11e') ' | ' ...
	strstop{1:end}]);
  if N < 102
     disp(['mean solution:' sprintf(' %+.1e', xmean)]);
     disp(['std deviation:' sprintf('  %.1e', sigma*sqrt(diagC))]);
     disp(sprintf('use plotcmaesdat.m for plotting the output at any time (option LogModulo must not be zero)'));
  end
  if exist('sfile', 'var') 
    disp(['Results saved in ' sfile]); 
  end
end

  out.arstopflags{irun} = stopflag;
  if any(strcmp(stopflag, 'fitness')) ...
	|| any(strcmp(stopflag, 'maxfunevals')) ...
	|| any(strcmp(stopflag, 'stoptoresume')) ...
	|| any(strcmp(stopflag, 'manual'))
  %  break; 
    break_result = 1;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE
newstate = state;
newstate.fmin = fmin;   newstate.xmin = xmin;
newstate.out = out;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE
  