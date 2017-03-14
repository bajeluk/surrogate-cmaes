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

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Algorithm comparison: CMA-ES, MA-ES, lmm-CMA-ES, saACMES, S-CMA-ES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimension 5.

data = {cmaes_data, ...
        maes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        scmaes_gp_data, ...
        scmaes_rf_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'MA-ES', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'S-CMA-ES GP', 'S-CMA-ES RF', 'DTS-CMA-ES'};

colors = [cmaesCol; maesCol; lmmCol; saacmesCol; scmaes_gpCol; scmaes_rfCol; dtsCol]/255;

plotFuns = 1:24;
plotDims = 5;

clear pdfNames
pdfNames = {};
for f = plotFuns
  for d = plotDims
    pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('alg_f%d_%dD', f, d));
  end
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', false, ...
                              'LineSpecification', {'-', '-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% Algorithm comparison: CMA-ES, MA-ES, lmm-CMA-ES, saACMES, S-CMA-ES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimension 20.

data = {cmaes_data, ...
        maes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        scmaes_gp_data, ...
        scmaes_rf_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'MA-ES', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'S-CMA-ES GP', 'S-CMA-ES RF', 'DTS-CMA-ES'};

colors = [cmaesCol; maesCol; lmmCol; saacmesCol; scmaes_gpCol; scmaes_rfCol; dtsCol]/255;

plotFuns = 1:24;
plotDims = 20;

clear pdfNames
pdfNames = {};
for f = plotFuns
  for d = plotDims
    pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('alg_f%d_%dD', f, d));
  end
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', false, ...
                              'LineSpecification', {'-', '-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% Aggregated algorithm comparison: CMA-ES, MA-ES, lmm-CMA-ES, saACMES, S-CMA-ES, DTS-CMA-ES  
% Aggregated  scaled function values in dimensions 5 and 20.

data = {cmaes_data, ...
        maes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        scmaes_gp_data, ...
        scmaes_rf_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'MA-ES', 'lmm-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'S-CMA-ES GP', 'S-CMA-ES RF', 'DTS-CMA-ES'};

colors = [cmaesCol; maesCol; lmmCol; saacmesCol; scmaes_gpCol; scmaes_rfCol; dtsCol]/255;

plotFuns = 1:24;
plotDims = [5, 20];

clear pdfNames
pdfNames = {};
for d = plotDims
  pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('alg_%dD', d));
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpecification', {'-', '-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)
                         
%% final clearing
close all
