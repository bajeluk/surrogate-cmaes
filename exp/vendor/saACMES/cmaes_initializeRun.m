function newstate = cmaes_initializeRuns(state)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE
[fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval,countevalNaN,irun ...
    flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals,stopMaxIter,stopFunEvals,stopIter, ... 
    stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues, ... 
    arrEqualFunvals, flgDiagonalOnly, flgActiveCMA,noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo, ... 
    flgscience,flgsaving,strsaving,flgsavingfinal, savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B,BD, ... 
    C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1,ccovmu,ccov1_sep, ccovmu_sep, ... 
    damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries,noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg, ... 
    stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop] = cmaes_loadState(state);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE


% ------------------------ Initialization -------------------------------

% Handle resuming of old run
flgresume = myevalbool(opts.Resume);
xmean = myeval(xstart); 
if all(size(xmean) > 1)
   xmean = mean(xmean, 2); % in case if xstart is a population
elseif size(xmean, 2) > 1
  xmean = xmean';
end 
if ~flgresume % not resuming a former run
  % Assign settings from input parameters and options for myeval...
  N = size(xmean, 1); numberofvariables = N; 
  lambda0 = floor(myeval(opts.PopSize) * myeval(opts.IncPopSize)^(irun-1)); 
  % lambda0 = floor(myeval(opts.PopSize) * 3^floor((irun-1)/2)); 
  popsize = lambda0;
  lambda = lambda0;
  insigma = myeval(insigma);
  if all(size(insigma) == [N 2]) 
    insigma = 0.5 * (insigma(:,2) - insigma(:,1));
  end
else % flgresume is true, do resume former run
  tmp = whos('-file', opts.SaveFilename);
  for i = 1:length(tmp)
    if strcmp(tmp(i).name, 'localopts');
      error('Saved variables include variable "localopts", please remove');
    end
  end
  local.opts = opts; % keep stopping and display options
  local.varargin = varargin;
  load(opts.SaveFilename); 
  varargin = local.varargin;
  flgresume = 1;

  % Overwrite old stopping and display options
  opts.StopFitness = local.opts.StopFitness; 
  %%opts.MaxFunEvals = local.opts.MaxFunEvals;
  %%opts.MaxIter = local.opts.MaxIter; 
  opts.StopFunEvals = local.opts.StopFunEvals; 
  opts.StopIter = local.opts.StopIter;  
  opts.TolX = local.opts.TolX;
  opts.TolUpX = local.opts.TolUpX;
  opts.TolFun = local.opts.TolFun;
  opts.TolHistFun = local.opts.TolHistFun;
  opts.StopOnStagnation = local.opts.StopOnStagnation; 
  opts.StopOnWarnings = local.opts.StopOnWarnings; 
  opts.ReadSignals = local.opts.ReadSignals; 
  opts.DispFinal = local.opts.DispFinal;
  opts.LogPlot = local.opts.LogPlot;
  opts.DispModulo = local.opts.DispModulo;
  opts.SaveVariables = local.opts.SaveVariables;
  opts.LogModulo = local.opts.LogModulo;
  opts.LogTime = local.opts.LogTime;
  clear local; % otherwise local would be overwritten during load
end
  
%--------------------------------------------------------------
% Evaluate options
stopFitness = myeval(opts.StopFitness); 
stopMaxFunEvals = myeval(opts.MaxFunEvals);  
stopMaxIter = myeval(opts.MaxIter);  
stopFunEvals = myeval(opts.StopFunEvals);  
stopIter = myeval(opts.StopIter);  
stopTolX = myeval(opts.TolX);
stopTolUpX = myeval(opts.TolUpX);
stopTolFun = myeval(opts.TolFun);
stopTolHistFun = myeval(opts.TolHistFun);
stopOnStagnation = myevalbool(opts.StopOnStagnation); 
stopOnWarnings = myevalbool(opts.StopOnWarnings); 
flgreadsignals = myevalbool(opts.ReadSignals);
flgWarnOnEqualFunctionValues = myevalbool(opts.WarnOnEqualFunctionValues);
flgEvalParallel = myevalbool(opts.EvalParallel);
stopOnEqualFunctionValues = myeval(opts.StopOnEqualFunctionValues);
arrEqualFunvals = zeros(1,10+N);
flgDiagonalOnly = myeval(opts.DiagonalOnly); 
flgActiveCMA = myeval(opts.CMA.active); 
noiseHandling = myevalbool(opts.Noise.on);
noiseMinMaxEvals = myeval(opts.Noise.minmaxevals);
noiseAlphaEvals = myeval(opts.Noise.alphaevals);
noiseCallback = myeval(opts.Noise.callback); 
flgdisplay = myevalbool(opts.DispFinal);
flgplotting = myevalbool(opts.LogPlot);
verbosemodulo = myeval(opts.DispModulo);
flgscience = myevalbool(opts.Science);
flgsaving = [];
strsaving = [];
if strfind(opts.SaveVariables, '-v6') 
  i = strfind(opts.SaveVariables, '%');
  if isempty(i) || i == 0 || strfind(opts.SaveVariables, '-v6') < i
    strsaving = '-v6';
    flgsaving = 1;
    flgsavingfinal = 1;
  end
end
if strncmp('final', opts.SaveVariables, 5)
  flgsaving = 0;
  flgsavingfinal = 1;
end
if isempty(flgsaving)
  flgsaving = myevalbool(opts.SaveVariables);
  flgsavingfinal = flgsaving;
end
savemodulo = myeval(opts.LogModulo);
savetime = myeval(opts.LogTime);

i = strfind(opts.LogFilenamePrefix, ' '); % remove everything after white space
if ~isempty(i)
  opts.LogFilenamePrefix = opts.LogFilenamePrefix(1:i(1)-1);
end

% TODO here silent option? set disp, save and log options to 0 

%--------------------------------------------------------------

if (isfinite(stopFunEvals) || isfinite(stopIter)) && ~flgsaving
  warning('To resume later the saving option needs to be set');
end


% Do more checking and initialization 
if flgresume % resume is on
  time.t0 = clock;
  if flgdisplay
    disp(['  resumed from ' opts.SaveFilename ]); 
  end
  if counteval >= stopMaxFunEvals 
    error(['MaxFunEvals exceeded, use StopFunEvals as stopping ' ...
	  'criterion before resume']);
  end
  if countiter >= stopMaxIter 
    error(['MaxIter exceeded, use StopIter as stopping criterion ' ...
	  'before resume']);
  end
  
else % flgresume
  % xmean = mean(myeval(xstart), 2); % evaluate xstart again, because of irun
  maxdx = myeval(opts.DiffMaxChange); % maximal sensible variable change
  mindx = myeval(opts.DiffMinChange); % minimal sensible variable change 
				      % can both also be defined as Nx1 vectors
  lbounds = myeval(opts.LBounds);		     
  ubounds = myeval(opts.UBounds);
  if length(lbounds) == 1
    lbounds = repmat(lbounds, N, 1);
  end
  if length(ubounds) == 1
    ubounds = repmat(ubounds, N, 1);
  end
  if isempty(insigma) % last chance to set insigma
    if all(lbounds > -Inf) && all(ubounds < Inf)
      if any(lbounds>=ubounds)
	error('upper bound must be greater than lower bound');
      end
      insigma = 0.3*(ubounds-lbounds);
      stopTolX = myeval(opts.TolX);  % reevaluate these
      stopTolUpX = myeval(opts.TolUpX);
    else
      error(['Initial step sizes (SIGMA) not determined']);
    end
  end

  % Check all vector sizes
  if size(xmean, 2) > 1 || size(xmean,1) ~= N
    error(['intial search point should be a column vector of size ' ...
	   num2str(N)]);
  elseif ~(all(size(insigma) == [1 1]) || all(size(insigma) == [N 1]))
    error(['input parameter SIGMA should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolX, 2) > 1 || ~ismember(size(stopTolX, 1), [1 N])
    error(['option TolX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolUpX, 2) > 1 || ~ismember(size(stopTolUpX, 1), [1 N])
    error(['option TolUpX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(maxdx, 2) > 1 || ~ismember(size(maxdx, 1), [1 N])
    error(['option DiffMaxChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(mindx, 2) > 1 || ~ismember(size(mindx, 1), [1 N])
    error(['option DiffMinChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(lbounds, 2) > 1 || ~ismember(size(lbounds, 1), [1 N])
    error(['option lbounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(ubounds, 2) > 1 || ~ismember(size(ubounds, 1), [1 N])
    error(['option ubounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  end
  
  % Initialize dynamic internal state parameters
  if any(insigma <= 0) 
    error(['Initial search volume (SIGMA) must be greater than zero']);
  end
  if max(insigma)/min(insigma) > 1e6
    error(['Initial search volume (SIGMA) badly conditioned']);
  end
  sigma = max(insigma);              % overall standard deviation
  pc = zeros(N,1); ps = zeros(N,1);  % evolution paths for C and sigma

  if length(insigma) == 1
    insigma = insigma * ones(N,1) ;
  end
  diagD = insigma/max(insigma);      % diagonal matrix D defines the scaling
  diagC = diagD.^2; 
  if flgDiagonalOnly ~= 1            % use at some point full covariance matrix
    B = eye(N,N);                      % B defines the coordinate system
    BD = B.*repmat(diagD',N,1);        % B*D for speed up only
    C = diag(diagC);                   % covariance matrix == BD*(BD)'
  end
  if flgDiagonalOnly
    B = 1; 
  end

  fitness.hist=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histsel=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histbest=[]; % history of fitness values
  fitness.histmedian=[]; % history of fitness values

  % Initialize boundary handling
  bnd.isactive = any(lbounds > -Inf) || any(ubounds < Inf); 
  if bnd.isactive
    if any(lbounds>ubounds)
      error('lower bound found to be greater than upper bound');
    end
    [xmean ti] = xintobounds(xmean, lbounds, ubounds); % just in case
    if any(ti)
      warning('Initial point was out of bounds, corrected');
    end
    bnd.weights = zeros(N,1);         % weights for bound penalty
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero 
    if bnd.flgscale ~= 0 
      bnd.scale = diagC/mean(diagC);
    else
      bnd.scale = ones(N,1);
    end
    
    idx = (lbounds > -Inf) | (ubounds < Inf);
    if length(idx) == 1
      idx = idx * ones(N,1);
    end
    bnd.isbounded = zeros(N,1);
    bnd.isbounded(find(idx)) = 1; 
    maxdx = min(maxdx, (ubounds - lbounds)/2);
    if any(sigma*sqrt(diagC) > maxdx)
      fac = min(maxdx ./ sqrt(diagC))/sigma;
      sigma = min(maxdx ./ sqrt(diagC));
      warning(['Initial SIGMA multiplied by the factor ' num2str(fac) ...
	       ', because it was larger than half' ...
	       ' of one of the boundary intervals']);
    end
    idx = (lbounds > -Inf) & (ubounds < Inf);
    dd = diagC;
    if any(5*sigma*sqrt(dd(idx)) < ubounds(idx) - lbounds(idx))
      warning(['Initial SIGMA is, in at least one coordinate, ' ...
	       'much smaller than the '...
	       'given boundary intervals. For reasonable ' ...
	       'global search performance SIGMA should be ' ...
	       'between 0.2 and 0.5 of the bounded interval in ' ...
	       'each coordinate. If all coordinates have ' ... 
	       'lower and upper bounds SIGMA can be empty']);
    end
    bnd.dfithist = 1;              % delta fit for setting weights
    bnd.aridxpoints = [];          % remember complete outside points
    bnd.arfitness = [];            % and their fitness
    bnd.validfitval = 0;
    bnd.iniphase = 1;
  end

  % ooo initial feval, for output only
  if irun == 1 
    out.solutions.bestever.x = xmean;
    out.solutions.bestever.f = Inf;  % for simpler comparison below
    out.solutions.bestever.evals = counteval;
    bestever = out.solutions.bestever;
  end
  if myevalbool(opts.EvalInitialX)
    fitness.hist(1)=feval(fitfun, xmean, varargin{:}); 
    fitness.histsel(1)=fitness.hist(1);
    counteval = counteval + 1;
    if fitness.hist(1) < out.solutions.bestever.f 
	out.solutions.bestever.x = xmean;
	out.solutions.bestever.f = fitness.hist(1);
	out.solutions.bestever.evals = counteval;
	bestever = out.solutions.bestever;
    end
  else
    fitness.hist(1)=NaN; 
    fitness.histsel(1)=NaN; 
  end
    
  % initialize random number generator
  if ischar(opts.Seed)
    randn('state', eval(opts.Seed));     % random number generator state
  else
    randn('state', opts.Seed);
  end
  %qqq
%  load(opts.SaveFilename, 'startseed');
%  randn('state', startseed);
%  disp(['SEED RELOADED FROM ' opts.SaveFilename]);
  startseed = randn('state');         % for retrieving in saved variables

  % Initialize further constants
  chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of 
				      %   ||N(0,I)|| == norm(randn(N,1))
  
  countiter = 0;
  % Initialize records and output
  if irun == 1
    time.t0 = clock;
    
    % TODO: keep also median solution? 
    out.evals = counteval;  % should be first entry
    out.stopflag = {};

    outiter = 0;

    % Write headers to output data files 
    filenameprefix = opts.LogFilenamePrefix; 
    if savemodulo && savetime
      filenames = {};
      filenames(end+1) = {'axlen'};
      filenames(end+1) = {'fit'};
      filenames(end+1) = {'stddev'};
      filenames(end+1) = {'xmean'};
      filenames(end+1) = {'xrecentbest'};
      str = [' (startseed=' num2str(startseed(2)) ...
             ', ' num2str(clock, '%d/%02d/%d %d:%d:%2.2f') ')'];
      for namecell = filenames(:)'
        name = namecell{:};

	[fid, err] = fopen(['./' filenameprefix name '.dat'], 'w');
	if fid < 1 % err ~= 0 
	  warning(['could not open ' filenameprefix name '.dat']);
	  filenames(find(strcmp(filenames,name))) = [];
	else
%	  fprintf(fid, '%s\n', ...
%	      ['<CMAES-OUTPUT version="' cmaVersion '">']);
%	  fprintf(fid, ['  <NAME>' name '</NAME>\n']);
%	  fprintf(fid, ['  <DATE>' date() '</DATE>\n']);
%	  fprintf(fid, '  <PARAMETERS>\n');
%	  fprintf(fid, ['    dimension=' num2str(N) '\n']);
%	  fprintf(fid, '  </PARAMETERS>\n');
	  % different cases for DATA columns annotations here
%	  fprintf(fid, '  <DATA');
	  if strcmp(name, 'axlen')
	     fprintf(fid, ['%%  columns="iteration, evaluation, sigma, ' ...
		 'max axis length, min axis length, ' ...
		 'all principal axes lengths (sorted square roots ' ...
                  'of eigenvalues of C)"' str]);
	  elseif strcmp(name, 'fit')
	    fprintf(fid, ['%%  columns="iteration, evaluation, sigma, axis ratio, bestever,' ...
		' best, median, worst fitness function value,' ...
		' further objective values of best"' str]);
	  elseif strcmp(name, 'stddev')
	    fprintf(fid, ['%%  columns=["iteration, evaluation, sigma, void, void, ' ...
		'stds==sigma*sqrt(diag(C))"' str]);
	  elseif strcmp(name, 'xmean')
	    fprintf(fid, ['%%  columns="iteration, evaluation, void, ' ...
                          'void, void, xmean"' str]);
	  elseif strcmp(name, 'xrecentbest')
	    fprintf(fid, ['%%  columns="iteration, evaluation, fitness, ' ...
                          'void, void, xrecentbest"' str]);
	  end
	  fprintf(fid, '\n'); % DATA
	  if strcmp(name, 'xmean')
	    fprintf(fid, '%ld %ld 0 0 0 ', 0, counteval); 
	    % fprintf(fid, '%ld %ld 0 0 %e ', countiter, counteval, fmean); 
%qqq	    fprintf(fid, msprintf('%e ', genophenotransform(out.genopheno, xmean)) + '\n'); 
	    fprintf(fid, '%e ', xmean);
            fprintf(fid, '\n'); 
	  end
	  fclose(fid); 
          clear fid; % preventing 
	end
      end % for files
    end % savemodulo
  end % irun == 1
  
end % else flgresume 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE
newstate = cmaes_saveState(state, fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval, ...
 countevalNaN,irun, flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals, ... 
 stopMaxIter,stopFunEvals,stopIter, stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,... 
 flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues,arrEqualFunvals, flgDiagonalOnly, flgActiveCMA, ... 
 noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo,flgscience,flgsaving, ... 
 strsaving,flgsavingfinal,savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B, ... 
 BD,C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1, ... 
 ccovmu,ccov1_sep, ccovmu_sep,damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries, ...
 noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg, stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE




