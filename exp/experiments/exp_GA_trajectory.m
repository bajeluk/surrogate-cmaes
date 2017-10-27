% EXAMPLE 15:  Continuous EDAs that learn mixtures of distributions
%              for  the trajectory problem (see previous examples for
%              details on this problem)

global MGADSMproblem;
thisDir = pwd();
cd('exp/vendor/mga'); % The spacecraft-trajectory problem instance directory
load('EdEdJ.mat');
cd(thisDir);
clear thisDir;

% The best found solution for 150.000 evaluations (can be improved...):
% x_min = [8400.66105421427, 6.85750764462776, 0.999342172061714, 0.50002169250583, 1114.15890187025, 2000, 0.514821614469194, 0.699088021681813, 3.79071058834825, 29.960723535454, -1.57102261983356, 1.06943051732632];
% f_min =279.509;
% worse one:
% x_min = [8510.07263718741, 5.12863867757723, 0.014467246215397, 0.501638295325626, 1421.88446602187, 1999.93550501781, 0.784333468763522, 0.794309742742636, 1.057696857128, 344.061516033378, 3.14136671748018, -2.76131610823915];
% f_min = 546.1515;

dim = 12;
PopSize = 5000; 
fitness = 'EvalSaga';
bounds(1,:) = [7000,0,0,0,50,300,0.01,0.01,1.05,8,-1*pi,-1*pi];
bounds(2,:) = [9100,7,1,1,2000,2000,0.90,0.90,7.00,500,pi,pi]; 
cache  = [0,0,1,0,1]; 

learning_params(1:5) = {'vars','ClusterPointsKmeans',10,'sqEuclidean',1};
edaparams{1} = {'learning_method','LearnMixtureofFullGaussianModels',learning_params};
edaparams{2} = {'sampling_method','SampleMixtureofFullGaussianModels',{PopSize,3}};
edaparams{3} = {'replacement_method','best_elitism',{'fitness_ordering'}};
selparams(1:2) = {0.1,'fitness_ordering'};
edaparams{4} = {'selection_method','truncation_selection',selparams};
edaparams{5} = {'repairing_method','SetWithinBounds_repairing',{}};
edaparams{6} = {'stop_cond_method','max_gen',{5000}};

maxfunevals = 150000;
ftarget = -Inf;
fitness = @(x) 1e6 - EvalSaga(x');
cmOptions = struct( ...
  'MaxFunEvals', min(1e8*dim, maxfunevals), ...
  'StopFitness', ftarget, ...
  'LBounds', bounds(1,:)', ...
  'UBounds', bounds(2,:)', ...
  'LogTime',  0, ...
  'SaveVariables', 'off', ...
  'LogModulo', 0, ...
  'Seed', 'inherit', ...
  'DispModulo', '10'); % , ...
  % 'PopSize', 20);

xstart = (bounds(1,:) + (0.1 + rand(1, dim) * 0.8) .* (bounds(2,:) - bounds(1,:)))';
sigmastart = ((bounds(2,:) - bounds(1,:)) / 2.1)';

% options = optimoptions('ga','PlotFcn',...
%     {@gaplotbestf,@gaplotbestindiv,@gaplotexpectation,@gaplotstopping});
options = optimoptions('ga', 'Display', 'iter');
[x_min, f_min, stopflag, output] = ga(fitness, dim, [], [], [], [], ...
    bounds(1,:), bounds(2,:), [], options);
