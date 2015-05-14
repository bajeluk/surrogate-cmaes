function state = cmaes( ...
    fitfun, ...    % name of objective/fitness function
    xstart, ...    % objective variables initial point, determines N
    insigma, ...   % initial coordinate wise standard deviation(s)
    inopts, ...    % options struct, see defopts below
    varargin )     % arguments passed to objective function 
% cmaes.m, Version 3.55.beta, last change: January, 2011 
% CMAES implements an Evolution Strategy with Covariance Matrix
% Adaptation (CMA-ES) for nonlinear function minimization.  For
% introductory comments and copyright (GPL) see end of file (type 
% 'type cmaes'). cmaes.m runs with MATLAB (Windows, Linux) and, 
% without data logging and plotting, it should run under Octave 
% (Linux, package octave-forge is needed).
%
% OPTS = CMAES returns default options. 
% OPTS = CMAES('defaults') returns default options quietly.
% OPTS = CMAES('displayoptions') displays options. 
% OPTS = CMAES('defaults', OPTS) supplements options OPTS with default 
% options. 
%
% XMIN = CMAES(FUN, X0, SIGMA[, OPTS]) locates the minimum XMIN of
% function FUN starting from column vector X0 with the initial
% coordinate wise search standard deviation SIGMA.
%
% Input arguments: 
%
%  FUN is a string function name like 'myfun'. FUN takes as argument a
%     column vector of size of X0 and returns a scalar. An easy way to
%     implement a hard non-linear constraint is to return NaN. Then,
%     this function evaluation is not counted and a newly sampled
%     point is tried immediately.
%
%   X0 is a column vector, or a matrix, or a string. If X0 is a matrix,
%     mean(X0, 2) is taken as initial point. If X0 is a string like
%     '2*rand(10,1)-1', the string is evaluated first.
%
%   SIGMA is a scalar, or a column vector of size(X0,1), or a string
%     that can be evaluated into one of these. SIGMA determines the
%     initial coordinate wise standard deviations for the search.
%     Setting SIGMA one third of the initial search region is
%     appropriate, e.g., the initial point in [0, 6]^10 and SIGMA=2
%     means cmaes('myfun', 3*rand(10,1), 2).  If SIGMA is missing and
%     size(X0,2) > 1, SIGMA is set to sqrt(var(X0')'). That is, X0 is
%     used as a sample for estimating initial mean and variance of the
%     search distribution.
%
%   OPTS (an optional argument) is a struct holding additional input
%     options. Valid field names and a short documentation can be
%     discovered by looking at the default options (type 'cmaes'
%     without arguments, see above). Empty or missing fields in OPTS
%     invoke the default value, i.e. OPTS needs not to have all valid
%     field names.  Capitalization does not matter and unambiguous
%     abbreviations can be used for the field names. If a string is
%     given where a numerical value is needed, the string is evaluated
%     by eval, where 'N' expands to the problem dimension
%     (==size(X0,1)) and 'popsize' to the population size. 
%
% [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = ...
%    CMAES(FITFUN, X0, SIGMA)
% returns the best (minimal) point XMIN (found in the last
% generation); function value FMIN of XMIN; the number of needed
% function evaluations COUNTEVAL; a STOPFLAG value as cell array,
% where possible entries are 'fitness', 'tolx', 'tolupx', 'tolfun',
% 'maxfunevals', 'maxiter', 'stoptoresume', 'manual',
% 'warnconditioncov', 'warnnoeffectcoord', 'warnnoeffectaxis',
% 'warnequalfunvals', 'warnequalfunvalhist', 'bug' (use
% e.g. any(strcmp(STOPFLAG, 'tolx')) or findstr(strcat(STOPFLAG,
% 'tolx')) for further processing); a record struct OUT with some
% more output, where the struct SOLUTIONS.BESTEVER contains the overall
% best evaluated point X with function value F evaluated at evaluation
% count EVALS. The last output argument BESTEVER equals 
% OUT.SOLUTIONS.BESTEVER. Moreover a history of solutions and 
% parameters is written to files according to the Log-options. 
%
% A regular manual stop can be achieved via the file signals.par. The
% program is terminated if the first two non-white sequences in any
% line of this file are 'stop' and the value of the LogFilenamePrefix
% option (by default 'outcmaes'). Also a run can be skipped.
% Given, for example, 'skip outcmaes run 2', skips the second run
% if option Restarts is at least 2, and another run will be started.
% 
% To run the code completely silently set Disp, Save, and Log options
% to 0.  With OPTS.LogModulo > 0 (1 by default) the most important
% data are written to ASCII files permitting to investigate the
% results (e.g. plot with function plotcmaesdat) even while CMAES is
% still running (which can be quite useful on expensive objective
% functions). When OPTS.SaveVariables==1 (default) everything is saved
% in file OPTS.SaveFilename (default 'variablescmaes.mat') allowing to
% resume the search afterwards by using the resume option.
%
% To find the best ever evaluated point load the variables typing
% "es=load('variablescmaes')" and investigate the variable
% es.out.solutions.bestever. 
%
% In case of a noisy objective function (uncertainties) set
% OPTS.Noise.on = 1. This option interferes presumably with some 
% termination criteria, because the step-size sigma will presumably
% not converge to zero anymore. If CMAES was provided with a 
% fifth argument (P1 in the below example, which is passed to the 
% objective function FUN), this argument is multiplied with the 
% factor given in option Noise.alphaevals, each time the detected 
% noise exceeds a threshold. This argument can be used within 
% FUN, for example, as averaging number to reduce the noise level.  
%
% OPTS.DiagonalOnly > 1 defines the number of initial iterations,
% where the covariance matrix remains diagonal and the algorithm has
% internally linear time complexity. OPTS.DiagonalOnly = 1 means
% keeping the covariance matrix always diagonal and this setting
% also exhibits linear space complexity. This can be particularly
% useful for dimension > 100. The default is OPTS.DiagonalOnly = 0. 
% 
% OPTS.CMA.active = 1 turns on "active CMA" with a negative update 
% of the covariance matrix and checks for positive definiteness. 
% OPTS.CMA.active = 2 does not check for pos. def. and is numerically
% faster. Active CMA usually speeds up the adaptation and might 
% become a default in near future. 
%
% The primary strategy parameter to play with is OPTS.PopSize, which
% can be increased from its default value.  Increasing the population
% size (by default linked to increasing parent number OPTS.ParentNumber)
% improves global search properties in exchange to speed. Speed
% decreases, as a rule, at most linearely with increasing population
% size. It is advisable to begin with the default small population
% size. The options Restarts and IncPopSize can be used for an
% automated multistart where the population size is increased by the
% factor IncPopSize (two by default) before each restart. X0 (given as
% string) is reevaluated for each restart. Stopping options
% StopFunEvals, StopIter, MaxFunEvals, and Fitness terminate the
% program, all others including MaxIter invoke another restart, where
% the iteration counter is reset to zero.
%
% Examples: 
%
%   XMIN = cmaes('myfun', 5*ones(10,1), 1.5); starts the search at
%   10D-point 5 and initially searches mainly between 5-3 and 5+3
%   (+- two standard deviations), but this is not a strict bound.
%   'myfun' is a name of a function that returns a scalar from a 10D
%   column vector.
%
%   opts.LBounds = 0; opts.UBounds = 10; 
%   X=cmaes('myfun', 10*rand(10,1), 5, opts);
%   search within lower bound of 0 and upper bound of 10. Bounds can
%   also be given as column vectors. If the optimum is not located
%   on the boundary, use rather a penalty approach to handle bounds. 
%
%   opts=cmaes; opts.StopFitness=1e-10;
%   X=cmaes('myfun', rand(5,1), 0.5, opts); stops the search, if
%   the function value is smaller than 1e-10.
%   
%   [X, F, E, STOP, OUT] = cmaes('myfun2', 'rand(5,1)', 1, [], P1, P2); 
%   passes two additional parameters to the function MYFUN2.
%
% See also FMINSEARCH, FMINUNC, FMINBND.


% TODO: 
%       write dispcmaesdat for Matlab (and Octave)
%       control savemodulo and plotmodulo via signals.par 


cmaVersion = '3.55.beta'; 

% ----------- Set Defaults for Input Parameters and Options -------------
% These defaults may be edited for convenience

% Input Defaults (obsolete, these are obligatory now)
definput.fitfun = 'felli'; % frosen; fcigar; see end of file for more
definput.xstart = rand(10,1); % 0.50*ones(10,1);
definput.sigma = 0.3;

% Options defaults: Stopping criteria % (value of stop flag)
defopts.StopFitness  = '-Inf % stop if f(xmin) < stopfitness, minimization';
defopts.MaxFunEvals  = 'Inf  % maximal number of fevals';
defopts.MaxIter      = '1e3*(N+5)^2/sqrt(popsize) % maximal number of iterations';
defopts.StopFunEvals = 'Inf  % stop after resp. evaluation, possibly resume later';
defopts.StopIter     = 'Inf  % stop after resp. iteration, possibly resume later';
defopts.TolX         = '1e-11*max(insigma) % stop if x-change smaller TolX';
defopts.TolUpX       = '1e3*max(insigma) % stop if x-changes larger TolUpX';
defopts.TolFun       = '1e-12 % stop if fun-changes smaller TolFun';
defopts.TolHistFun   = '1e-13 % stop if back fun-changes smaller TolHistFun';
defopts.StopOnStagnation = 'on  % stop when fitness stagnates for a long time';
defopts.StopOnWarnings = 'yes  % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.StopOnEqualFunctionValues = '2 + N/3  % number of iterations';  

% Options defaults: Other
defopts.DiffMaxChange = 'Inf  % maximal variable change(s), can be Nx1-vector';
defopts.DiffMinChange = '0    % minimal variable change(s), can be Nx1-vector';
defopts.WarnOnEqualFunctionValues = ...
    'yes  % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.LBounds = '-Inf % lower bounds, scalar or Nx1-vector'; 
defopts.UBounds = 'Inf  % upper bounds, scalar or Nx1-vector'; 
defopts.EvalParallel = 'no   % objective function FUN accepts NxM matrix, with M>1?';
defopts.EvalInitialX = 'yes  % evaluation of initial solution';
defopts.Restarts     = '0    % number of restarts ';
defopts.IncPopSize   = '2    % multiplier for population size before each restart';

defopts.PopSize      = '(4 + floor(3*log(N)))  % population size, lambda'; 
defopts.ParentNumber = 'floor(popsize/2)       % AKA mu, popsize equals lambda';
defopts.RecombinationWeights = 'superlinear decrease % or linear, or equal';
defopts.DiagonalOnly = '0*(1+100*N/sqrt(popsize))+(N>=1000)  % C is diagonal for given iterations, 1==always'; 
defopts.Noise.on = '0  % uncertainty handling is off by default'; 
defopts.Noise.reevals  = '1*ceil(0.05*lambda)  % nb. of re-evaluated for uncertainty measurement';
defopts.Noise.theta = '0.5  % threshold to invoke uncertainty treatment'; % smaller: more likely to diverge
defopts.Noise.cum = '0.3  % cumulation constant for uncertainty'; 
defopts.Noise.cutoff = '2*lambda/3  % rank change cutoff for summation';
defopts.Noise.alphasigma  = '1+2/(N+10)  % factor for increasing sigma'; % smaller: slower adaptation
defopts.Noise.epsilon  = '1e-7  % additional relative perturbation before reevaluation';
defopts.Noise.minmaxevals = '[1 inf]  % min and max value of 2nd arg to fitfun, start value is 5th arg to cmaes'; 
defopts.Noise.alphaevals  = '1+2/(N+10)  % factor for increasing 2nd arg to fitfun'; 
defopts.Noise.callback = '[]  % callback function when uncertainty threshold is exceeded';
% defopts.TPA = 0; 
defopts.CMA.cs = '(mueff+2)/(N+mueff+3)  % cumulation constant for step-size'; 
   %qqq cs = (mueff^0.5)/(N^0.5+mueff^0.5) % the short time horizon version
defopts.CMA.damps = '1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs  % damping for step-size';
defopts.CMA.ccum = '4/(N+4)  % cumulation constant for covariance matrix'; 
defopts.CMA.ccov1 = '2 / ((N+1.3)^2+mueff)  % learning rate for rank-one update'; 
defopts.CMA.ccovmu = '2 * (mueff-2+1/mueff) / ((N+2)^2+mueff) % learning rate for rank-mu update'; 
defopts.CMA.on     = 'yes'; 
defopts.CMA.active = '0  % active CMA 1: neg. updates with pos. def. check, 2: neg. updates'; 

flg_future_setting = 0;  % testing for possible future variant(s)
if flg_future_setting    
  disp('in the future')

  % damps setting from Brockhoff et al 2010
  %   this damps diverges with popsize 400:
  %   cmaeshtml('benchmarkszero', ones(20,1)*2, 5, o, 15);
  defopts.CMA.damps = '2*mueff/lambda + 0.3 + cs  % damping for step-size';  % cs: for large mueff
  % how about:
  % defopts.CMA.damps = '2*mueff/lambda + 0.3 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs  % damping for step-size';

  % ccum adjusted for large mueff, better on schefelmult? 
  % TODO: this should also depend on diagonal option!?  
  defopts.CMA.ccum = '(4 + mueff/N) / (N+4 + 2*mueff/N)  % cumulation constant for pc';

  defopts.CMA.active = '1  % active CMA 1: neg. updates with pos. def. check, 2: neg. updates'; 
end
  
defopts.Resume   = 'no   % resume former run from SaveFile'; 
defopts.Science  = 'on  % off==do some additional (minor) problem capturing, NOT IN USE'; 
defopts.ReadSignals = 'on  % from file signals.par for termination, yet a stumb';
defopts.Seed = 'sum(100*clock)  % evaluated if it is a string';
defopts.DispFinal  = 'on   % display messages like initial and final message';
defopts.DispModulo = '100  % [0:Inf], disp messages after every i-th iteration';
defopts.SaveVariables = 'on   % [on|final|off][-v6] save variables to .mat file';
defopts.SaveFilename = 'variablescmaes.mat  % save all variables, see SaveVariables'; 
defopts.LogModulo = '1    % [0:Inf] if >1 record data less frequently after gen=100';
defopts.LogTime   = '25   % [0:100] max. percentage of time for recording data';
defopts.LogFilenamePrefix = 'outcmaes  % files for output data'; 
defopts.LogPlot = 'off    % plot while running using output data files';

%qqqkkk 
%defopts.varopt1 = ''; % 'for temporary and hacking purposes'; 
%defopts.varopt2 = ''; % 'for temporary and hacking purposes'; 
defopts.UserData = 'for saving data/comments associated with the run';
defopts.UserDat2 = ''; 'for saving data/comments associated with the run';

% ---------------------- Handling Input Parameters ----------------------

if nargin < 1 || isequal(fitfun, 'defaults') % pass default options
  if nargin < 1
    disp('Default options returned (type "help cmaes" for help).');
  end
  xmin = defopts;
  if nargin > 1 % supplement second argument with default options
    xmin = getoptions(xstart, defopts);
  end
  state = xmin;
  return;
end

if isequal(fitfun, 'displayoptions')
 names = fieldnames(defopts); 
 for name = names'
   disp([name{:} repmat(' ', 1, 20-length(name{:})) ': ''' defopts.(name{:}) '''']); 
 end
 return; 
end

input.fitfun = fitfun; % record used input
if isempty(fitfun)
  % fitfun = definput.fitfun; 
  % warning(['Objective function not determined, ''' fitfun ''' used']);
  error(['Objective function not determined']);
end
if ~ischar(fitfun)
  error('first argument FUN must be a string');
end


if nargin < 2 
  xstart = [];
end

input.xstart = xstart;
if isempty(xstart)
  % xstart = definput.xstart;  % objective variables initial point
  % warning('Initial search point, and problem dimension, not determined');
  error('Initial search point, and problem dimension, not determined');
end

if nargin < 3 
  insigma = [];
end
if isa(insigma, 'struct')
  error(['Third argument SIGMA must be (or eval to) a scalar '...
	   'or a column vector of size(X0,1)']);
end
input.sigma = insigma;
if isempty(insigma)
  if all(size(myeval(xstart)) > 1)
    insigma = std(xstart, 0, 2); 
    if any(insigma == 0)
      error(['Initial search volume is zero, choose SIGMA or X0 appropriate']);
    end
  else
    % will be captured later
    % error(['Initial step sizes (SIGMA) not determined']);
  end
end

% Compose options opts
if nargin < 4 || isempty(inopts) % no input options available
  inopts = []; 
  opts = defopts;
else
  opts = getoptions(inopts, defopts);
end
i = strfind(opts.SaveFilename, ' '); % remove everything after white space
if ~isempty(i)
  opts.SaveFilename = opts.SaveFilename(1:i(1)-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
counteval = 0; countevalNaN = 0; 
irun = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE
state = cmaes_saveState({}, fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting, ... 
 nargin,input,opts,counteval, countevalNaN,irun, [], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[], ... 
 [],[],[],[],[],[],[],[],[], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[], ... 
 [],[],[],[],[],[],[],[],[], [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%
% Example:
%   In the source-code of the function that needs optional
%   arguments, the last argument is the struct of optional
%   arguments:
%
%   function results = myfunction(mandatory_arg, inopts)
%     % Define four default options
%     defopts.PopulationSize = 200;
%     defopts.ParentNumber = 50;
%     defopts.MaxIterations = 1e6;
%     defopts.MaxSigma = 1;
%  
%     % merge default options with input options
%     opts = getoptions(inopts, defopts);
%
%     % Thats it! From now on the values in opts can be used
%     for i = 1:opts.PopulationSize
%       % do whatever
%       if sigma > opts.MaxSigma
%         % do whatever
%       end
%     end
%   
%   For calling the function myfunction with default options:
%   myfunction(argument1, []);
%   For calling the function myfunction with modified options:
%   opt.pop = 100; % redefine PopulationSize option
%   opt.PAR = 10;  % redefine ParentNumber option
%   opt.maxiter = 2; % opt.max=2 is ambiguous and would result in an error 
%   myfunction(argument1, opt);

%
% 04/07/19: Entries can be structs itself leading to a recursive
%           call to getoptions. 
%

if nargin < 2 || isempty(defopts) % no default options available
  opts=inopts;
  return;
elseif isempty(inopts) % empty inopts invoke default options
  opts = defopts;
  return;
elseif ~isstruct(defopts) % handle a single option value
  if isempty(inopts) 
    opts = defopts;
  elseif ~isstruct(inopts)
    opts = inopts;
  else
    error('Input options are a struct, while default options are not');
  end
  return;
elseif ~isstruct(inopts) % no valid input options
  error('The options need to be a struct or empty');
end

  opts = defopts; % start from defopts 
  % if necessary overwrite opts fields by inopts values
  defnames = fieldnames(defopts);
  idxmatched = []; % indices of defopts that already matched
  for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
      for i = 1:size(defnames, 1)
	idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
      end
    else
	idx = strncmpi(defnames, name, length(name));
    end
    if sum(idx) > 1
      error(['option "' name '" is not an unambigous abbreviation. ' ...
	     'Use opts=RMFIELD(opts, ''' name, ...
	     ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
      defname  = defnames{find(idx)}; 
      if ismember(find(idx), idxmatched)
	error(['input options match more than ones with "' ...
	       defname '". ' ...
	       'Use opts=RMFIELD(opts, ''' name, ...
	       ''') to remove the field from the struct.']);
      end
      idxmatched = [idxmatched find(idx)];
      val = getfield(inopts, name);
      % next line can replace previous line from MATLAB version 6.5.0 on and in octave
      % val = inopts.(name);
      if isstruct(val) % valid syntax only from version 6.5.0
	opts = setfield(opts, defname, ...
	    getoptions(val, getfield(defopts, defname))); 
      elseif isstruct(getfield(defopts, defname)) 
      % next three lines can replace previous three lines from MATLAB 
      % version 6.5.0 on
      %   opts.(defname) = ...
      %      getoptions(val, defopts.(defname)); 
      % elseif isstruct(defopts.(defname)) 
	warning(['option "' name '" disregarded (must be struct)']); 
      elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
	opts = setfield(opts, defnames{find(idx)}, val);
	% next line can replace previous line from MATLAB version 6.5.0 on
	% opts.(defname) = inopts.(name); 
      end
    else
      warning(['option "' name '" disregarded (unknown field name)']);
    end
  end

% CHANGES
% 10/11/11: (3.52.beta) boundary handling: replace max with min in change 
%           rate formula. Active CMA: check of pos.def. improved. 
%           Plotting: value of lambda appears in the title. 
% 10/04/03: (3.51.beta) active CMA cleaned up. Equal fitness detection
%           looks into history now. 
% 10/03/08: (3.50.beta) "active CMA" revised and bug-fix of ambiguous
%           option Noise.alpha -> Noise.alphasigma. 
% 09/10/12: (3.40.beta) a slightly modified version of "active CMA", 
%           that is a negative covariance matrix update, use option
%           CMA.active. In 10;30;90-D the gain on ftablet is a factor 
%           of 1.6;2.5;4.4 (the scaling improves by sqrt(N)). On 
%           Rosenbrock the gain is about 25%. On sharp ridge the 
%           behavior is improved. Cigar is unchanged. 
% 09/08/10: local plotcmaesdat remains in backround  
% 09/08/10: bug-fix in time management for data writing, logtime was not 
%        considered properly (usually not at all). 
% 09/07/05: V3.24: stagnation termination added 
% 08/09/27: V3.23: momentum alignment is out-commented and de-preciated
% 08/09/25: V3.22: re-alignment of sigma and C was buggy
% 08/07/15: V3.20, CMA-parameters are options now. ccov and mucov were replaced
%        by ccov1 \approx ccov/mucov and ccovmu \approx (1-1/mucov)*ccov
% 08/06/30: file name xrecent was change to xrecentbest (compatible with other 
%        versions)
% 08/06/29: time stamp added to output files 
% 08/06/28: bug fixed with resume option, commentary did not work 
% 08/06/28: V3.10, uncertainty (noise) handling added (re-implemented), according
%        to reference "A Method for Handling Uncertainty..." from below.
% 08/06/28: bug fix: file xrecent was empty 
% 08/06/01: diagonalonly clean up. >1 means some iterations. 
% 08/05/05: output is written to file preventing an increasing data
%        array and ease long runs. 
% 08/03/27: DiagonalOnly<0 learns for -DiagonalOnly iterations only the 
%        diagonal with a larger learning rate. 
% 08/03 (2.60): option DiagonalOnly>=1 invokes a time- and space-linear  
%        variant with only diagonal elements of the covariance matrix 
%        updating.  This can be useful for large dimensions, say > 100. 
% 08/02: diag(weights) * ... replaced with repmat(weights,1,N) .* ...
%        in C update, implies O(mu*N^2) instead of O(mu^2*N + mu*N^2). 
% 07/09: tolhistfun as termination criterion added, "<" changed to
%        "<=" also for TolFun to allow for stopping on zero difference. 
%        Name tolfunhist clashes with option tolfun. 
% 07/07: hsig threshold made slighly smaller for large dimension, 
%        useful for lambda < lambda_default. 
% 07/06: boundary handling: scaling in the boundary handling
%        is omitted now, see bnd.flgscale. This seems not to
%        have a big impact. Using the scaling is worse on rotated
%        functions, but better on separable ones. 
% 07/05: boundary handling: weight i is not incremented anymore
%        if xmean(i) moves towards the feasible space. Increment
%        factor changed to 1.2 instead of 1.1. 
% 07/05: boundary handling code simplified not changing the algorithm
% 07/04: bug removed for saving in octave
% 06/11/10: more testing of outcome of eig, fixed max(D) to max(diag(D))
% 06/10/21: conclusive final bestever assignment in the end 
% 06/10/21: restart and incpopsize option implemented for restarts
%        with increasing population size, version 2.50. 
% 06/09/16: output argument bestever inserted again for convenience and
%        backward compatibility
% 06/08: output argument out and struct out reorganized. 
% 06/01: Possible parallel evaluation included as option EvalParallel
% 05/11: Compatibility to octave implemented, package octave-forge
%   is needed. 
% 05/09: Raise of figure and waiting for first plots improved
% 05/01: Function coordinatesystem cleaned up. 
% 05/01: Function prctile, which requires the statistics toolbox,
%        replaced by myprctile. 
% 05/01: Option warnonequalfunctionvalues included. 
% 04/12: Decrease of sigma removed. Problems on fsectorsphere can 
%        be addressed better by adding search space boundaries. 
% 04/12: Boundary handling simpyfied. 
% 04/12: Bug when stopping criteria tolx or tolupx are vectors. 
% 04/11: Three input parameters are obligatory now. 
% 04/11: Bug in boundary handling removed: Boundary weights can decrease now. 
% 04/11: Normalization for boundary weights scale changed. 
% 04/11: VerboseModulo option bug removed. Documentation improved. 
% 04/11: Condition for increasing boundary weights changed.
% 04/10: Decrease of sigma when fitness is getting consistenly
%        worse. Addresses the problems appearing on fsectorsphere for
%        large population size.
% 04/10: VerboseModulo option included. 
% 04/10: Bug for condition for increasing boundary weights removed.
% 04/07: tolx depends on initial sigma to achieve scale invariance
%        for this stopping criterion. 
% 04/06: Objective function value NaN is not counted as function
%        evaluation and invokes resampling of the search point. 
% 04/06: Error handling for eigenvalue beeing zero (never happens
%        with default parameter setting)
% 04/05: damps further tuned for large mueff 
%      o Details for stall of pc-adaptation added (variable hsig 
%        introduced). 
% 04/05: Bug in boundary handling removed: A large initial SIGMA was
%        corrected not until *after* the first iteration, which could
%        lead to a complete failure.
% 04/05: Call of function range (works with stats toolbox only) 
%        changed to myrange. 
% 04/04: Parameter cs depends on mueff now and damps \propto sqrt(mueff)
%        instead of \propto mueff. 
%      o Initial stall to adapt C (flginiphase) is removed and
%        adaptation of pc is stalled for large norm(ps) instead.
%      o Returned default options include documentation. 
%      o Resume part reorganized.
% 04/03: Stopflag becomes cell-array. 

% ---------------------------------------------------------------
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
% nonlinear function minimization. To be used under the terms of the
% GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
% Author (copyright): Nikolaus Hansen, 2001-2008. 
% e-mail: nikolaus.hansen AT inria.fr
% URL:http://www.bionik.tu-berlin.de/user/niko
% References: See below. 
% ---------------------------------------------------------------
%
% GENERAL PURPOSE: The CMA-ES (Evolution Strategy with Covariance
% Matrix Adaptation) is a robust search method which should be
% applied, if derivative based methods, e.g. quasi-Newton BFGS or
% conjucate gradient, (supposably) fail due to a rugged search
% landscape (e.g. noise, local optima, outlier, etc.). On smooth
% landscapes CMA-ES is roughly ten times slower than BFGS. For up to
% N=10 variables even the simplex direct search method (Nelder & Mead)
% is often faster, but far less robust than CMA-ES.  To see the
% advantage of the CMA, it will usually take at least 30*N and up to
% 300*N function evaluations, where N is the search problem dimension.
% On considerably hard problems the complete search (a single run) is
% expected to take at least 30*N^2 and up to 300*N^2 function
% evaluations.
%
% SOME MORE COMMENTS: 
% The adaptation of the covariance matrix (e.g. by the CMA) is
% equivalent to a general linear transformation of the problem
% coding. Nevertheless every problem specific knowlegde about the best
% linear transformation should be exploited before starting the
% search. That is, an appropriate a priori transformation should be
% applied to the problem. This also makes the identity matrix as
% initial covariance matrix the best choice.
%
% The strategy parameter lambda (population size, opts.PopSize) is the
% preferred strategy parameter to play with.  If results with the
% default strategy are not satisfactory, increase the population
% size. (Remark that the crucial parameter mu (opts.ParentNumber) is
% increased proportionally to lambda). This will improve the
% strategies capability of handling noise and local minima. We
% recomment successively increasing lambda by a factor of about three,
% starting with initial values between 5 and 20. Casually, population
% sizes even beyond 1000+100*N can be sensible.
%
%
% ---------------------------------------------------------------
%%% REFERENCES
%
% The equation numbers refer to 
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer. 
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% 
% Further references:
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%
% Ros, R. and N. Hansen (2008). A Simple Modification in CMA-ES
% Achieving Linear Time and Space Complexity. To appear in Tenth
% International Conference on Parallel Problem Solving from Nature
% PPSN X, Proceedings, Berlin: Springer.
%
% Hansen, N., A.S.P. Niederberger, L. Guzzella and P. Koumoutsakos
% (2009?). A Method for Handling Uncertainty in Evolutionary
% Optimization with an Application to Feedback Control of
% Combustion. To appear in IEEE Transactions on Evolutionary
% Computation.


