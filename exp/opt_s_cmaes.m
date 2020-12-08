function [x, y_evals, stopflag, varargout] = opt_s_cmaes(FUN, dim, ftarget, maxfunevals, id, varargin)
% minimizes FUN in DIM dimensions by multistarts of fminsearch.
% ftarget and maxfunevals are additional external termination conditions,
% where at most 2 * maxfunevals function evaluations are conducted.
% fminsearch was modified to take as input variable usual_delta to
% generate the first simplex.
% set options, make sure we always terminate
% with restarts up to 2*maxfunevals are allowed
%
% last varargin argument can be 'exppath' -- path to the experiment's data

% Be aware: 'id' is an additional parameter!

fDelta = 1e-8;

cmOptions = struct( ...
  'MaxFunEvals', min(1e8*dim, maxfunevals), ...
  'StopFitness', ftarget, ...
  'LBounds', -5, ...
  'UBounds',  5, ...
  'LogTime',  0, ...
  'SaveVariables', 'off', ...
  'LogModulo', 0, ...
  'Seed', 'inherit', ...
  'DispModulo', '10');

y_evals = [];

if (nargin >= 6)  exppath = [varargin{1} filesep];
  else              exppath = sgParams.experimentPath;  end

load([exppath 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');
[bbParams, sgParams, cmParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

for fname = fieldnames(cmParams)'
  cmOptions.(fname{1}) = cmParams.(fname{1});
end

if (nargin >= 7)  xstart = varargin{2};
  else              xstart = 8 * rand(dim, 1) - 4;  end
if (nargin >= 8)  datapath = varargin{3};
  else              datapath = exppath;  end
if (nargin >= 9)  iinstance = varargin{4};
  else              iinstance = NaN;  end

  % Info about tested function is for debugging purposes
  bbob_handlesF = benchmarks('handles');
  noisyHandles = benchmarksnoisy('handles');
  extraHandles = benchmarksextra('handles');
  bbob_handlesF(100+(1:length(noisyHandles))) = noisyHandles;
  bbob_handlesF(200+(1:length(extraHandles))) = extraHandles;
  sgParams.modelOpts.bbob_func = bbob_handlesF{bbParams.functions(1)};
  sgParams.expFileID = [num2str(bbParams.functions(1)) '_' num2str(dim) 'D_' num2str(id)];
  sgParams.instance  = iinstance;
  sgParams.fopt      = ftarget - fDelta;
  sgParams.datapath  = datapath;
  [~, sgParams.exp_id] = fileparts(sgParams.experimentPath);
  % DEBUG: generate data for testing model regresssion
  % TODO: comment this line! :)
  % sgParams.saveModelTrainingData = [ 10 25 50 100 200 300 470 700 900 1200 1500 2000 2400 ];

  [x, fmin, counteval, stopflag, out, bestever, y_eval] = s_cmaes(FUN, xstart, 8/3, cmOptions, 'SurrogateOptions', sgParams);

  n_y_evals = size(y_eval,1);
  y_eval(:,1) = y_eval(:,1) - (ftarget - fDelta) * ones(n_y_evals,1);
  y_evals = [y_evals; y_eval];

  if (nargout > 3)
    varargout{1} = out;
  else
    varargout = cell(0);
  end

end % function
