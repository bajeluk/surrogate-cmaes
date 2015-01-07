
exp_id = 'exp_pure_cmaes_01';
exp_description = 'Pure CMA-ES';

machines = {'u-pl28'};
login = 'bajel3am';
matlabcommand = 'matlab_mff_2014b'; % '/afs/ms/@sys/bin/matlab';
logMatlabOutput = true;

% BBOB parameters
bbParamDef(1).name   = 'dimensions';
bbParamDef(1).values = {2, 5, 10};
bbParamDef(2).name   = 'functions';
bbParamDef(2).values = {2, 3, 8};
% dimensions  = [10];     % which dimensions to optimize, subset of [2 3 5 10 20 40];
% functions   = [8];      % function ID's to optimize (2 Sphere, 3 Rastrigin, 8 Rosenbrock)
bbParamDef(3).name   = 'opt_function';
bbParamDef(3).values = {@opt_s_cmaes};
% opt_function = @opt_s_cmaes;    % function being optimized -- BBOB wrap-around with header
%                                 % xbest = function( fun, dim, ftarget, maxfunevals )
bbParamDef(4).name   = 'instances';
bbParamDef(4).values = {[1:3]}; % {[1:5 31:40]};  % default is [1:5, 31:40]
bbParamDef(5).name   = 'maxfunevals';   % MAXFUNEVALS - 10*dim is a short test-experiment
bbParamDef(5).values = {'250 * dim'};   % increment maxfunevals successively
                                
% Surrogate model parameter lists
sgParamDef(1).name   = 'evoControl';
sgParamDef(1).values = {'none'};

% CMA-ES parameters
cmParamDef(1).name   = 'PopSize';
cmParamDef(1).values = {'(4 + floor(3*log(N)))'}; %, '(8 + floor(6*log(N)))'};
cmParamDef(2).name   = 'Restarts';
cmParamDef(2).values = {4};

% path to current file -- do not change this
pathstr = fileparts(mfilename('fullpath'));
exppath  = [pathstr filesep exp_id];
exppath_short  = pathstr;
[s,mess,messid] = mkdir(exppath);
[s,mess,messid] = mkdir([exppath filesep 'cmaes_results']);
addpath(exppath);
save([exppath filesep 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');

% run the rest of the scripts generation
generateShellScripts
