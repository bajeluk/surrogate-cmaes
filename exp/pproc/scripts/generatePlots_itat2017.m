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

adts_he50_path      = fullfile(exppath, 'exp_doubleEC_23_adapt05_v3');
adts_le50_he75_path = fullfile(exppath, 'exp_doubleEC_23_adapt05_v2');
adts_le75_he75_path = fullfile(exppath, 'exp_doubleEC_23_adapt05');

dts_path = fullfile(exppath, 'exp_doubleEC_23');
cmaes_path = fullfile(exppath, 'CMA-ES');
saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
lmm_path = fullfile(exppath, 'lmm-CMA-ES');

% load data
dataFolders = {adts_he50_path; ...
               adts_le50_he75_path; ...
               adts_le75_he75_path; ...
               dts_path; ...
               cmaes_path; ...
               saacmes_path; ...
               lmm_path};
             
[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings

% adaptive DTS settings id
clear findSet
% quantiles: lowErr 0.5, highErr 0.5
findSet.DTAdaptive_lowErr = 0.5;
findSet.DTAdaptive_highErr = 0.5;
findSet.DTAdaptive_updateRate = 0.3;
adts_le50_he50_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le50_he50_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le50_he50_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.5, highErr 0.75
findSet.DTAdaptive_lowErr = 0.5;
findSet.DTAdaptive_highErr = 0.75;
findSet.DTAdaptive_updateRate = 0.3;
adts_le50_he75_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le50_he75_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le50_he75_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.75, highErr 0.5
findSet.DTAdaptive_lowErr = 0.75;
findSet.DTAdaptive_highErr = 0.5;
findSet.DTAdaptive_updateRate = 0.3;
adts_le75_he50_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le75_he50_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le75_he50_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.75, highErr 0.75
findSet.DTAdaptive_lowErr = 0.75;
findSet.DTAdaptive_highErr = 0.75;
findSet.DTAdaptive_updateRate = 0.3;
adts_le75_he75_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le75_he75_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le75_he75_ur5_Id = getStructIndex(settings, findSet);

% DTS settings id
clear findSet
findSet.DTAdaptive_defaultErr = [];
dts_Id = getStructIndex(settings, findSet);

% CMA-ES, saACMES, lmm-CMA-ES settings id
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

adts_le50_he50_ur3_data = evals(:, :, adts_le50_he50_ur3_Id);
adts_le50_he50_ur4_data = evals(:, :, adts_le50_he50_ur4_Id);
adts_le50_he50_ur5_data = evals(:, :, adts_le50_he50_ur5_Id);

adts_le50_he75_ur3_data = evals(:, :, adts_le50_he75_ur3_Id);
adts_le50_he75_ur4_data = evals(:, :, adts_le50_he75_ur4_Id);
adts_le50_he75_ur5_data = evals(:, :, adts_le50_he75_ur5_Id);

adts_le75_he50_ur3_data = evals(:, :, adts_le75_he50_ur3_Id);
adts_le75_he50_ur4_data = evals(:, :, adts_le75_he50_ur4_Id);
adts_le75_he50_ur5_data = evals(:, :, adts_le75_he50_ur5_Id);

adts_le75_he75_ur3_data = evals(:, :, adts_le75_he75_ur3_Id);
adts_le75_he75_ur4_data = evals(:, :, adts_le75_he75_ur4_Id);
adts_le75_he75_ur5_data = evals(:, :, adts_le75_he75_ur5_Id);

% color settings
cmaesCol     = getAlgColors('cmaes');
saacmesCol   = getAlgColors('saacmes');
dtsCol       = getAlgColors('dtscmaes');
lmmCol       = getAlgColors('lmmcmaes');

adts_le50_he50_ur3_col = getAlgColors(1);
adts_le50_he50_ur4_col = getAlgColors(2);
adts_le50_he50_ur5_col = getAlgColors(3);

adts_le50_he75_ur3_col = getAlgColors(4);
adts_le50_he75_ur4_col = getAlgColors(5);
adts_le50_he75_ur5_col = getAlgColors(6);

adts_le75_he50_ur3_col = getAlgColors(7);
adts_le75_he50_ur4_col = getAlgColors(8);
adts_le75_he50_ur5_col = getAlgColors(9);

adts_le75_he75_ur3_col = getAlgColors(10);
adts_le75_he75_ur4_col = getAlgColors(11);
adts_le75_he75_ur5_col = getAlgColors(12);

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

%% ADTS-CMA-ES settings comparison: lowErr, highErr, updateRate
% Aggregation of function values across dimensions 2, 3, 5, 10, 20.

close all

set_data = {adts_le50_he50_ur3_data, ...
            adts_le50_he50_ur4_data, ...
            adts_le50_he50_ur5_data, ...
            adts_le50_he75_ur3_data, ...
            adts_le50_he75_ur4_data, ...
            adts_le50_he75_ur5_data, ...
            adts_le75_he50_ur3_data, ...
            adts_le75_he50_ur4_data, ...
            adts_le75_he50_ur5_data, ...
            adts_le75_he75_ur3_data, ...
            adts_le75_he75_ur4_data, ...
            adts_le75_he75_ur5_data};

set_datanames = {'$\\err_\\text{min}$, $\\err_\\text{max}$'};

tableFunc = funcSet.BBfunc;
tableDims = funcSet.dims;

resultTable = fullfile(tableFolder, 'rankTable.tex');
      
[table, ranks] = rankingTable(set_data, 'DataNames', set_datanames, ...
                           'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                           'TableFuns', tableFunc, 'TableDims', tableDims,...
                           'Evaluations', [20 40 80], ...
                           'ResultFile', resultTable);

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
