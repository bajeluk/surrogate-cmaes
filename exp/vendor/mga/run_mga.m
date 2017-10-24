% For a description of the spacecraft trajectory problem and links
% to the problem instances see (Vinko_et_al:2007) 

global MGADSMproblem
load 'EdEdJ.mat'  % This instance can be downloaded from the ESA web page 
NumbVar = 12;
PopSize = 1000; 
F = 'sagas'; % The function should be modified to enforce maximization, i.e g(x) = -1*f(x)
% F = 'EvalSaga'; % The function should be modified to enforce maximization, i.e g(x) = -1*f(x)
Card(1,:) = [7000,0,0,0,50,300,0.01,0.01,1.05,8,-1*pi,-1*pi];
Card(2,:) = [9100,7,1,1,2000,2000,0.90,0.90,7.00,500,pi,pi]; 
cache  = [0,0,0,0,0];  
learning_params(1:5) = {'vars','ClusterPointsKmeans',10,'sqEuclidean',1};
edaparams{1} = {'learning_method','LearnMixtureofFullGaussianModels',learning_params};
edaparams{2} = {'sampling_method','SampleMixtureofFullGaussianModels',{PopSize,3}};
edaparams{3} = {'replacement_method','best_elitism',{'fitness_ordering'}};
selparams(1:2) = {0.1,'fitness_ordering'};
edaparams{4} = {'selection_method','truncation_selection',selparams};
edaparams{5} = {'repairing_method','SetWithinBounds_repairing',{}};
edaparams{6} = {'stop_cond_method','max_gen',{5000}};

% [AllStat,Cache]=RunEDA(PopSize,NumbVar,F,Card,cache,edaparams) 

