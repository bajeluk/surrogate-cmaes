% runs the timing experiment for MY_OPTIMIZER. fgeneric.m
% and benchmarks.m must be in the path of MATLAB/Octave

exp_id = 'exp_geneEC_08_10D';
dimensions = [2 5 10 20];
id = 1; % ... for CMA-ES options
maxfunevals = '100*dim';

pathstr = fileparts(mfilename('fullpath'));
exppath = [pathstr filesep '..' filesep '..' filesep 'experiments' filesep exp_id];

load([exppath filesep 'scmaes_params.mat']);
bbParams = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
MY_OPTIMIZER = bbParams.opt_function;

addpath(pathstr);  % should point to fgeneric.m etc.

more off;  % in octave pagination is on by default

timings = [];
runs = [];
dims = [];
for dim = dimensions
  nbrun = 0;
  ftarget = fgeneric('initialize', 8, 1, 'tmp');
  tic;
  while toc < 400  % at least 30 seconds
    MY_OPTIMIZER(@fgeneric, dim, ftarget, myeval(maxfunevals), id);  % adjust maxfunevals
    nbrun = nbrun + 1;
  end  % while
  timings(end+1) = toc / fgeneric('evaluations');
  dims(end+1) = dim;    % not really needed
  runs(end+1) = nbrun;  % not really needed
  fgeneric('finalize');
  disp([['Dimensions:' sprintf(' %11d ', dims)]; ...
        ['      runs:' sprintf(' %11d ', runs)]; ...
        [' times [s]:' sprintf(' %11.1e ', timings)]]);
end

