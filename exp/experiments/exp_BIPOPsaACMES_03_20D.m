exp_id = 'exp_BIPOPsaACMES_03_20D';
exp_description = 'Loshchilov''s saACM-ES algorithm (non-BIPOP), 24 functions, 15 instances, 20D, 3rd generation experiment';

machines = {'machine1'};

login = 'bajel3am';
if (strfind(mfilename('fullpath'), 'afs'))
  matlabcommand = '/afs/ms/@sys/bin/matlab';
else
  matlabcommand = 'matlab_mff_2014b';
end
logMatlabOutput = true;
logDir = '/storage/plzen1/home/bajeluk/public';

% BBOB parameters
bbParamDef(1).name   = 'dimensions';
bbParamDef(1).values = {20};      % {2, 5 10};
bbParamDef(2).name   = 'functions';
bbParamDef(2).values = num2cell(1:24);  % {1, 2, 3, 5, 6, 8, 10, 11, 12, 13, 14, 20, 21};
% dimensions  = [10];     % which dimensions to optimize, subset of [2 3 5 10 20 40];
% functions   = [8];      % function ID's to optimize (2 Sphere, 3 Rastrigin, 8 Rosenbrock)
bbParamDef(3).name   = 'opt_function';
bbParamDef(3).values = {@opt_xacmes};
% opt_function = @opt_s_cmaes;    % function being optimized -- BBOB wrap-around with header
%                                 % xbest = function( fun, dim, ftarget, maxfunevals )
bbParamDef(4).name   = 'instances';
bbParamDef(4).values = {[1:5 31:40]}; % 31:40]};   % default is [1:5, 31:40]
bbParamDef(5).name   = 'maxfunevals';   % MAXFUNEVALS - 10*dim is a short test-experiment
bbParamDef(5).values = {'250 * dim'};   % increment maxfunevals successively
                                
% Surrogate model parameter lists
sgParamDef(1).name   = 'lambdaMult';
lambdaMult = 1;
% lambdaMult = ones(1,40);
% lambdaMult([2,3,5]) = 1;
% lambdaMult([10]) = 10;
% lambdaMult([20]) = 100;
% lambdaMult([40]) = 1000;
sgParamDef(1).values = { lambdaMult };
% BIPOP = false && useSCMAES = false && withSurr = 1 
% stands for Ilya's BIPOP-saACM-ES
sgParamDef(2).name   = 'BIPOP';
sgParamDef(2).values = { 1 };
sgParamDef(3).name   = 'useSCMAES';
sgParamDef(3).values = { 0 };
sgParamDef(4).name   = 'withSurr';
sgParamDef(4).values = { 1 };

% CMA-ES parameters
cmParamDef(1).name   = 'PopSize';
cmParamDef(1).values = {'(4 + floor(3*log(N)))'}; %, '(8 + floor(6*log(N)))'};
cmParamDef(2).name   = 'Restarts';
cmParamDef(2).values = {100};

% path to current file -- do not change this
pathstr = fileparts(mfilename('fullpath'));
exppath  = [pathstr filesep exp_id];
exppath_short  = pathstr;
[s,mess,messid] = mkdir(exppath);
[s,mess,messid] = mkdir([exppath filesep 'cmaes_results']);
addpath(exppath);

% Save the directory of the experiment data for debugging purposes
sgParamDef(end+1).name = 'experimentPath';
sgParamDef(end).values = { exppath };

save([exppath filesep 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef', 'exp_id', 'exppath_short', 'logDir');

% run the rest of the scripts generation
generateShellScriptsMetacentrum
