%% GECCO 2017 workshop plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values of compared algorithms.
% 
% Created for GECCO 2017 workshop article.

%% load data

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
% TODO: replace by more recent results
dts_path = fullfile(exppath, 'DTS-CMA-ES_05_2pop');
maes_path = fullfile(exppath, 'MA-ES');

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
findSet.algName = 'CMA-ES';
cma_Id = getStructIndex(settings, findSet);
findSet.algName = 'MA-ES';
ma_Id = getStructIndex(settings, findSet);
findSet.algName = 'BIPOP-saACM-k';
saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'DTS-CMA-ES_05_2pop';
dts_Id = getStructIndex(settings, findSet);
findSet.algName = 'lmm-CMA-ES';
lmm_Id = getStructIndex(settings, findSet);
             
% extract data
scmaes_gp_data = evals(:, :, gp_Id);
scmaes_rf_data = evals(:, :, rf_Id);
cmaes_data     = evals(:, :, cma_Id);
maes_data      = evals(:, :, ma_Id);
saacmes_data   = evals(:, :, saacm_Id);
dtscmaes_data  = evals(:, :, dts_Id);
lmmcmaes_data  = evals(:, :, lmm_Id);

% color settings
scmaes_rfCol = getAlgColors(1);
scmaes_gpCol = getAlgColors('scmaes');
maesCol      = getAlgColors(2);
cmaesCol     = getAlgColors('cmaes');
saacmesCol   = getAlgColors('saacmes');
dtsCol       = getAlgColors('dtscmaes');
lmmCol       = getAlgColors('lmmcmaes');

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Algorithm comparison: DTS-CMA-ES, S-CMA-ES, saACMES, SMAC, CMA-ES
% Aggregation of function values across dimensions 2, 5, 10, 20.

data = {cmaes_data, ...
        maes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        scmaes_gp_data, ...
        scmaes_rf_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'MA-ES', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'S-CMA-ES GP', 'S-CMA-ES RF', 'DTS-CMA-ES'};

colors = [cmaesCol; maesCol; lmmCol; saacmesCol; scmaes_gpCol; scmaes_rfCol; dtsCol]/255;

plotDims = [2, 5, 10, 20];

clear pdfNames
pdfNames = fullfile(plotResultsFolder, 'alg2_5_10_20D');

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', true, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpecification', {'-.', '-.', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', 100);

                              
                            
% print2pdf(han, pdfNames, 1)
                         
%% final clearing
close all
