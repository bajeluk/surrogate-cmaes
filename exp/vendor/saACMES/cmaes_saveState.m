function state = cmaes_saveState(state, fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval, ...
 countevalNaN,irun, flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals, ... 
 stopMaxIter,stopFunEvals,stopIter, stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,... 
 flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues,arrEqualFunvals, flgDiagonalOnly, flgActiveCMA, ... 
 noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo,flgscience,flgsaving, ... 
 strsaving,flgsavingfinal,savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B, ... 
 BD,C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1, ... 
 ccovmu,ccov1_sep, ccovmu_sep,damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries, ...
 noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg,stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop)


state.fitfun = fitfun;
state.xstart = xstart;
state.insigma = insigma;
state.inopts = inopts;
state.varargin = varargin;
state.cmaVersion = cmaVersion;
state.definput = definput;
state.defopts = defopts;
state.flg_future_setting = flg_future_setting;
state.nargin = nargin;
state.input = input;
state.opts = opts;
state.counteval = counteval;
state.countevalNaN = countevalNaN;
state.irun = irun;
state.flgresume = flgresume;
state.xmean = xmean;
state.N = N;
state.numberofvariables = numberofvariables;
state.lambda0 = lambda0;
state.popsize = popsize;
state.lambda = lambda;
state.lambda_last = lambda_last;
state.stopFitness = stopFitness;
state.stopMaxFunEvals = stopMaxFunEvals;
state.stopMaxIter = stopMaxIter;
state.stopFunEvals = stopFunEvals;
state.stopIter = stopIter;
state.stopTolFun = stopTolFun;
state.stopTolHistFun = stopTolHistFun;
state.stopOnStagnation = stopOnStagnation;
state.stopOnWarnings = stopOnWarnings;
state.flgreadsignals = flgreadsignals;
state.flgWarnOnEqualFunctionValues = flgWarnOnEqualFunctionValues;
state.flgEvalParallel = flgEvalParallel;
state.stopOnEqualFunctionValues = stopOnEqualFunctionValues;
state.arrEqualFunvals = arrEqualFunvals;
state.flgDiagonalOnly = flgDiagonalOnly;
state.flgActiveCMA = flgActiveCMA;
state.noiseHandling = noiseHandling;
state.noiseMinMaxEvals = noiseMinMaxEvals;
state.noiseAlphaEvals = noiseAlphaEvals;
state.noiseCallback = noiseCallback;
state.flgdisplay = flgdisplay;
state.flgplotting = flgplotting;
state.verbosemodulo = verbosemodulo;
state.flgscience = flgscience;
state.flgsaving = flgsaving;
state.strsaving = strsaving;
state.flgsavingfinal = flgsavingfinal;
state.savemodulo = savemodulo;
state.savetime = savetime;
state.time = time;
state.maxdx = maxdx;
state.mindx = mindx;
state.lbounds = lbounds;
state.ubounds = ubounds;
state.stopTolX = stopTolX;
state.stopTolUpX = stopTolUpX;
state.sigma = sigma;
state.pc = pc;
state.diagD = diagD;
state.diagC = diagC;
state.B = B;
state.BD = BD;
state.C = C;
state.fitness = fitness;
state.bnd = bnd;
state.out = out;
state.startseed = startseed;
state.chiN = chiN;
state.countiter = countiter;
state.outiter = outiter;
state.filenameprefix = filenameprefix;
state.filenames = filenames;
state.lambda_hist = lambda_hist;
state.mu = mu;
state.weights = weights;
state.mueff = mueff;
state.cc = cc;
state.cs = cs;
state.ccov1 = ccov1;
state.ccovmu = ccovmu;
state.ccov1_sep = ccov1_sep;
state.ccovmu_sep = ccovmu_sep;
state.damps = damps;
state.noiseReevals = noiseReevals; 
state.noiseAlpha = noiseAlpha;
state.noiseEpsilon = noiseEpsilon;
state.noiseTheta = noiseTheta;
state.noisecum = noisecum;
state.noiseCutOff = noiseCutOff;
state.arx = arx;
state.arxvalid = arxvalid;
state.tries = tries;
state.noiseS = noiseS;
state.noiseSS = noiseSS;
state.noiseN = noiseN;
state.xold = xold;
state.zmean = zmean;
state.fmean = fmean;
state.ps = ps;
state.neg = neg;
state.stopflag = stopflag;
state.noiseX = noiseX;
state.iterplotted = iterplotted;

state.arfitness = arfitness;
state.Xnew_sorted = Xnew_sorted;
state.invsqrtC = invsqrtC;

state.stop = stop;