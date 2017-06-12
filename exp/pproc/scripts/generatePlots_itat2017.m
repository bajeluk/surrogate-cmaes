%% ITAT 2017 plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values of compared algorithms.
% 
% Created for GECCO 2017 workshop article.

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'itat2017_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 5, 10, 20];
maxEvals = 250;
  
% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'itat2017paper');
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

adts_path = fullfile(exppath, 'DTS-CMA-ES_05_2pop');
dts_path = fullfile(exppath, 'exp_doubleEC_23');

cmaes_path = fullfile(exppath, 'CMA-ES');
saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
lmm_path = fullfile(exppath, 'lmm-CMA-ES');

% load data
dataFolders = {adts_path; ...
               dts_path; ...
               cmaes_path; ...
               saacmes_path; ...
               lmm_path};
             
[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings

% adaptive DTS settings id

clear findSet
findSet.evoControlRestrictedParam = 0.05;
findSet.modelOpts.predictionType = 'sd2';
dts_Id = getStructIndex(settings, findSet);

clear findSet
findSet.algName = 'CMA-ES';
cma_Id = getStructIndex(settings, findSet);
findSet.algName = 'BIPOP-saACM-k';
saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'lmm-CMA-ES';
lmm_Id = getStructIndex(settings, findSet);
             
% extract data
cmaes_data     = evals(:, :, cma_Id);
saacmes_data   = evals(:, :, saacm_Id);
dtscmaes_data  = evals(:, :, dts_Id);
lmmcmaes_data  = evals(:, :, lmm_Id);

% color settings
cmaesCol     = getAlgColors('cmaes');
saacmesCol   = getAlgColors('saacmes');
dtsCol       = getAlgColors('dtscmaes');
lmmCol       = getAlgColors('lmmcmaes');

% aggregate data & settings
data = {cmaes_data, ...
        lmmcmaes_data, ...
        saacmes_data, ...
        dtscmaes_data};

datanames = {'CMA-ES', 'lmm-CMA-ES', '{}^{s*}ACMES-k', 'DTS-CMA-ES'};

colors = [cmaesCol; lmmCol; saacmesCol; dtsCol]/255;

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimension 5.

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
                              'LegendOption', 'first', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% Algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimension 20.

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
                              'LegendOption', 'first', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% Aggregated algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Aggregated  scaled function values in dimensions 5 and 20.

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

%% Multiple comparison of algorithms with a statistical posthoc test.

close all

tableFunc = funcSet.BBfunc;
tableDims = [5, 20];

resultDuelTable = fullfile(tableFolder, 'duelTable.tex');

datanames = {'CMA-ES', 'lmm-CMA-ES', '\\saACMES-k', 'DTS-CMA-ES'};

[table, ranks] = duelTable(data, 'DataNames', datanames, ...
                            'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                            'TableFuns', tableFunc, 'TableDims', tableDims, ...
                            'Evaluations', [1/3 1], ...
                            'ResultFile', resultDuelTable);

%% final clearing
close all
