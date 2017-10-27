function [x, y_evals, stopflag, varargout] = opt_EGNA(FUN, dim, ftarget, maxfunevals, id, varargin)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed
%
% varargin/1 argument can be 'xstart' -- starting point for the optimization
% varargin/2 can be 'cmOptions' -- options for the CEDA algorithm


  fDelta = 1e-8;

  cmOptions = struct( ...
    'LBounds', -5, ...
    'UBounds',  5, ...
    'PopSize', '200*dim', ...
    'ClusterPointsKmeans', 10, ...
    'sqEuclidean', 1, ...
    'SelectionRatio', 0.1);

  y_evals = [];

  if (nargin >= 6), exppath = [varargin{1} filesep];
    else              exppath = '';  end
  if (nargin >= 7), xstart = varargin{2};
    else              xstart = 8 * rand(dim, 1) - 4;  end

  load([exppath 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');
  [~, ~, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

  for fname = fieldnames(cmParams)'
    cmOptions.(fname{1}) = cmParams.(fname{1});
  end

  lb = cmOptions.LBounds * ones(dim, 1);
  ub = cmOptions.UBounds * ones(dim, 1);
  bounds = [lb'; ub'];
  if (isstr(FUN))
    FUN = str2func(FUN);
  end
  fitness = @(x) -FUN(x');

  cache  = [0,1,0,0,1];  
  learning_params = {'vars','ClusterPointsKmeans',cmOptions.ClusterPointsKmeans, ...
    'sqEuclidean',cmOptions.sqEuclidean};
  edaparams{1} = {'learning_method','LearnMixtureofFullGaussianModels',learning_params};
  edaparams{2} = {'sampling_method','SampleMixtureofFullGaussianModels',{myeval(cmOptions.PopSize),3}};
  edaparams{3} = {'replacement_method','best_elitism',{'fitness_ordering'}};
  selparams(1:2) = {cmOptions.SelectionRatio,'fitness_ordering'};
  edaparams{4} = {'selection_method','truncation_selection',selparams};
  edaparams{5} = {'repairing_method','SetWithinBounds_repairing',{}};
  edaparams{6} = {'stop_cond_method','maxgen_maxval', ...
      {ceil(maxfunevals/myeval(cmOptions.PopSize)), ftarget}};
  edaparams{7} = {'verbose_method', 'none', {}};

  % EGNA itself
  %
  [AllStat,Cache]=RunEDA(myeval(cmOptions.PopSize),dim,fitness,bounds,cache,edaparams);

  fmins = cell2num(AllStat(:,1));
  fmins = fmins(1:5:end);
  evals = cell2num(AllStat(:,5));
  
  x = AllStat{end, 2};
  y_evals = [fmins, evals, NaN(length(fmins),3)];

  stopflag = [];
  if (nargout > 3)
    varargout{1} = [];
  else
    varargout = cell(0);
  end
end % function
