% Startup script for Surrogate CMA-ES

disp('executing surrogate-cmaes startup script...');

addpath(genpath(fullfile(pwd, 'src')));
addpath(genpath(fullfile(pwd, 'test')));
addpath(fullfile(pwd, 'exp'));
addpath(genpath(fullfile(pwd, 'exp/pproc')));
addpath(genpath(fullfile(pwd, 'exp/util')));
addpath(genpath(fullfile(pwd, 'exp/vendor')));
addpath(genpath(fullfile(pwd, 'exp/log')));
addpath(fullfile(pwd, 'exp/experiments'));

% run('src/vendor/gpml-matlab-v3.2/startup.m');
run('src/vendor/gpml_v4.0/startup.m');
