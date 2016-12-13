function [x, y_evals, stopflag, varargout] = opt_cmaes(FUN, DIM, ftarget, maxfunevals, id, varargin)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed
%
% last varargin argument can be 'exppath' -- path to the experiment's data


fDelta = 1e-8;

cmOptions = struct( ...
  'MaxFunEvals', min(1e8*DIM, maxfunevals), ...
  'StopFitness', ftarget, ...
  'LBounds', -5, ...
  'UBounds',  5, ...
  'LogTime',  0, ...
  'SaveVariables', 'off', ...
  'Seed', 'inherit', ...
  'LogModulo', 0, ...
  'DispModulo', '10');

y_evals = [];

if (nargin >= 6)  exppath = [varargin{1} filesep];
  else              exppath = '';  end
if (nargin >= 7)  xstart = varargin{2};
  else              xstart = 8 * rand(dim, 1) - 4;  end

load([exppath 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');
[~, ~, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

for fname = fieldnames(cmParams)'
  cmOptions.(fname{1}) = cmParams.(fname{1});
end
  cmOptions.PopSize = '(4 + floor(3*log(N)))';  % CMA-ES default

  [x, fmin, counteval, stopflag, out, bestever, y_eval] = s_cmaes(FUN, xstart, 8/3, cmOptions);

  n_y_evals = size(y_eval,1);
  y_eval(:,1) = y_eval(:,1) - (ftarget - fDelta) * ones(n_y_evals,1);
  y_evals = [y_evals; y_eval];

  if (nargout > 3)
    varargout{1} = out;
  else
    varargout = cell(0);
  end

end % function
