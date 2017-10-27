% EXAMPLE 15:  Continuous EDAs that learn mixtures of distributions
%              for  the trajectory problem (see previous examples for
%              details on this problem)

global MGADSMproblem;
thisDir = pwd();
cd('exp/vendor/mga'); % The spacecraft-trajectory problem instance directory
load('EdEdJ.mat');
cd(thisDir);
clear thisDir;

dim = 12;
PopSize = 5000; 
fitness = 'EvalSaga';
bounds(1,:) = [7000,0,0,0,50,300,0.01,0.01,1.05,8,-1*pi,-1*pi];
bounds(2,:) = [9100,7,1,1,2000,2000,0.90,0.90,7.00,500,pi,pi]; 
cache  = [0,0,1,0,1]; 

maxfunevals = 150000;
ftarget = 0;
fitness = @(x) -1e6 + EvalSaga(x);

cmOptions = struct( ...
  'PopSize', '200*dim', ...
  'CopulaType', 'Gaussian', ...
  'SelectionRatio', 0.1);

learning_params = {cmOptions.CopulaType};
edaparams{1} = {'learning_method','LearnCopulaModel',learning_params};
edaparams{2} = {'sampling_method','SampleCopulaModel',{myeval(cmOptions.PopSize),3}};
edaparams{3} = {'replacement_method','best_elitism',{'fitness_ordering'}};
selparams(1:2) = {cmOptions.SelectionRatio,'fitness_ordering'};
edaparams{4} = {'selection_method','truncation_selection',selparams};
edaparams{5} = {'repairing_method','SetWithinBounds_repairing',{}};
edaparams{6} = {'stop_cond_method','maxgen_maxval', ...
    {ceil(maxfunevals/myeval(cmOptions.PopSize)), ftarget}};
% edaparams{7} = {'verbose_method', 'none', {}};

[AllStat,Cache] = RunEDA(myeval(cmOptions.PopSize),dim,fitness,bounds,cache,edaparams);

fmins = cell2num(AllStat(:,1));
fmins = fmins(1:5:end);
evals = cell2num(AllStat(:,5));

x_min = AllStat{end, 2};
f_min = fmins(end);
y_evals = [fmins, evals, NaN(length(fmins),3)];
