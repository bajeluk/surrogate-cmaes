%% GECCO 2017 workshop article plots
% Graphs show the dependence of minimal function
% values on the number of function values of compared algorithms.
%
% Let $\Delta_f$ be the minimal distance found from the function optimum
% for the considered number of fitness function evaluations.
% The graphs depict a scaled logarithm of $\Delta_f$
% depending on FE/D. Since all the algorithms ran for each function and dimension
% on 15 independent instances, only the empirical medians~$\Delta_f^\text{med}$ over those 15 runs of $\Delta_f$ were taken
% for further processing.
% The scaled logarithms of $\Delta_f^\textrm{med}$ are calculated as
% 
% $$ \Delta_f^{\log} = \frac{\log \Delta_f^\textrm{med} -
% \Delta_f^\textrm{MIN}}{\Delta_f^\textrm{MAX} - \Delta_f^\textrm{MIN}} 
% \log_{10} \left(1 / 10^{-8}\right) + \log_{10} 10^{-8} \,,$$
%
% where $\Delta_f^\textrm{MIN}$
% ($\Delta_f^\textrm{MAX}$) is the minimum (maximum) 
% $\log \Delta_f^\textrm{med}$ found among all the compared algorithms
% for the particular function $f$ and dimension $D$ between 0 and 
% 250 FE/D. Afterwards, graphs of $\Delta_f^{\log}$ can be aggregated across
% arbitrary number of functions and dimensions. 
% 
% Created for GECCO 2017 workshop article.

%%

% Load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'gecco2017workshop_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 5, 10, 20];
maxEvals = 250;
  
% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2017workshop');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
if ~isdir(plotResultsFolder)
  mkdir(plotResultsFolder)
end
if ~isdir(tableFolder)
  mkdir(tableFolder)
end

% path settings
exppath = fullfile('exp', 'experiments');

gen_path = fullfile(exppath, 'exp_geneEC_10');
gen_path20D = fullfile(exppath, 'exp_geneEC_10_20D');
dts_path = fullfile(exppath, 'DTS-CMA-ES_05_2pop');
maes_path = fullfile(exppath, 'exp_maesEC_14_2_10_cmaes_20D');

cmaes_path = fullfile(exppath, 'CMA-ES');
saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
lmm_path = fullfile(exppath, 'lmm-CMA-ES');

% load data
dataFolders = {gen_path; ...
               gen_path20D; ...
               dts_path; ...
               maes_path; ...
               cmaes_path; ...
               saacmes_path; ...
               lmm_path};
             
[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings
clear findSet
findSet.modelType = 'gp';
findSet.evoControlModelGenerations = 5;
gp_Id = getStructIndex(settings, findSet);

findSet.modelType = 'rf';
findSet.evoControlModelGenerations = 1;
rf_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'fvalues';
ma_Id = getStructIndex(settings, findSet);

clear findSet
findSet.algName = 'CMA-ES';
cma_Id = getStructIndex(settings, findSet);
findSet.algName = 'BIPOP-saACM-k';
saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'lmm-CMA-ES';
lmm_Id = getStructIndex(settings, findSet);
findSet.algName = 'DTS-CMA-ES_05_2pop';
dts_Id = getStructIndex(settings, findSet);
             
% extract data
scmaes_gp_data = evals(:, :, gp_Id);
scmaes_rf_data = evals(:, :, rf_Id);
cmaes_data     = evals(:, :, cma_Id);
maes_data      = evals(:, :, ma_Id);
saacmes_data   = evals(:, :, saacm_Id);
dtscmaes_data  = evals(:, :, dts_Id);
lmmcmaes_data  = evals(:, :, lmm_Id);

% color settings
scmaes_rfCol = getAlgColors(2);
scmaes_gpCol = getAlgColors('scmaes');
maesCol      = getAlgColors(3);
cmaesCol     = getAlgColors('cmaes');
saacmesCol   = getAlgColors('saacmes');
dtsCol       = getAlgColors('dtscmaes');
lmmCol       = getAlgColors('lmmcmaes');

% aggregate data & settings
data = {cmaes_data, ...
        maes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        scmaes_gp_data, ...
        scmaes_rf_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'MA-ES', 'lmm-CMA-ES', '{}^{s*}ACMES-k', 'S-CMA-ES GP', 'S-CMA-ES RF', 'DTS-CMA-ES'};

colors = [cmaesCol; maesCol; lmmCol; saacmesCol; scmaes_gpCol; scmaes_rfCol; dtsCol]/255;

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Algorithm comparison: CMA-ES, MA-ES, lmm-CMA-ES, saACMES, S-CMA-ES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimensions 2, 5, 10, and 20.

for f = funcSet.BBfunc

  %% 
  close all

  fprintf('Function %d\n', f)
  han = relativeFValuesPlot(data, ...
                              'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, ...
                              'PlotFuns', f, ...
                              'AggregateDims', false, ...
                              'AggregateFuns', true, ...
                              'DataNames', datanames, ...
                              'Colors', colors, ...
                              'FunctionNames', true, ...
                              'LegendOption', 'out', ...
                              'MaxEvals', maxEvals, ...
                              'Statistic', @median);
end

%% Summary graphs
% Summary graphs are averaged through all functions for individual algorithms in separate dimensionalities.

close all

han = relativeFValuesPlot(data, ...
                            'DataDims', funcSet.dims, ...
                            'DataFuns', funcSet.BBfunc, ...
                            'PlotFuns', funcSet.BBfunc, ...
                            'AggregateDims', false, ...
                            'AggregateFuns', true, ...
                            'DataNames', datanames, ...
                            'Colors', colors, ...
                            'FunctionNames', true, ...
                            'LegendOption', 'out', ...
                            'MaxEvals', maxEvals, ...
                            'Statistic', @median);

%%

% Final clearing
close all
clear s f