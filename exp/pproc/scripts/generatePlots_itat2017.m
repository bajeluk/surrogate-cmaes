%% ITAT 2017 plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values of compared algorithms.
% 
% Created for GECCO 2017 workshop article.

%% load data

fprintf('Loading data...\n')

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'itat2017_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10, 20];
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
findSet.DTAdaptive_updateRateDown = 'obj.updateRate';
findSet.DTAdaptive_maxRatio = 1;

% quantiles: lowErr 0.5, highErr 0.5
findSet.DTAdaptive_lowErr = '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.11; -0.0092; -0.13; 0.044; 0.14]';
findSet.DTAdaptive_highErr = '@(x) [ones(size(x,1),1) x(:,1) x(:,2) x(:,1).*x(:,2) x(:,2).^2] * [0.18; -0.0027; 0.44; 0.0032; -0.14]';
findSet.DTAdaptive_updateRate = 0.3;
adts_le50_he50_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le50_he50_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le50_he50_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.5, highErr 0.75
findSet.DTAdaptive_lowErr = '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.11; -0.0092; -0.13; 0.044; 0.14]';
findSet.DTAdaptive_highErr = '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.35; -0.047; 0.44; 0.044; -0.19]';
findSet.DTAdaptive_updateRate = 0.3;
adts_le50_he75_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le50_he75_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le50_he75_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.75, highErr 0.5
findSet.DTAdaptive_lowErr = '@(x) [ones(size(x,1),1) x(:,1) x(:,2) x(:,1).*x(:,2) x(:,2).^2] * [0.17; -0.00067; -0.095; 0.0087; 0.15]';
findSet.DTAdaptive_highErr = '@(x) [ones(size(x,1),1) x(:,1) x(:,2) x(:,1).*x(:,2) x(:,2).^2] * [0.18; -0.0027; 0.44; 0.0032; -0.14]';
findSet.DTAdaptive_updateRate = 0.3;
adts_le75_he50_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le75_he50_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le75_he50_ur5_Id = getStructIndex(settings, findSet);

% quantiles: lowErr 0.75, highErr 0.75
findSet.DTAdaptive_lowErr = '@(x) [ones(size(x,1),1) x(:,1) x(:,2) x(:,1).*x(:,2) x(:,2).^2] * [0.17; -0.00067; -0.095; 0.0087; 0.15]';
findSet.DTAdaptive_highErr = '@(x) [ones(size(x,1),1) log(x(:,1)) x(:,2) log(x(:,1)).*x(:,2) x(:,2).^2] * [0.35; -0.047; 0.44; 0.044; -0.19]';
findSet.DTAdaptive_updateRate = 0.3;
adts_le75_he75_ur3_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.4;
adts_le75_he75_ur4_Id = getStructIndex(settings, findSet);

findSet.DTAdaptive_updateRate = 0.5;
adts_le75_he75_ur5_Id = getStructIndex(settings, findSet);

% DTS settings id
clear findSet
findSet.modelOpts.inputFraction = 1;
findSet.modelOpts.predictionType = 'sd2';
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

adts_le50_he50_ur3_col = getAlgColors(3);
adts_le50_he50_ur4_col = getAlgColors(2);
adts_le50_he50_ur5_col = getAlgColors(1);

adts_le50_he75_ur3_col = getAlgColors(6);
adts_le50_he75_ur4_col = getAlgColors(5);
adts_le50_he75_ur5_col = getAlgColors(4);

adts_le75_he50_ur3_col = getAlgColors(9);
adts_le75_he50_ur4_col = getAlgColors(8);
adts_le75_he50_ur5_col = getAlgColors(7);

adts_le75_he75_ur3_col = getAlgColors(12);
adts_le75_he75_ur4_col = getAlgColors(11);
adts_le75_he75_ur5_col = getAlgColors(10);

% aggregate data & settings
data = {...
  adts_le50_he50_ur3_data, ...
  adts_le50_he75_ur3_data, ...
  adts_le75_he50_ur3_data, ...
  adts_le75_he75_ur3_data, ...
  cmaes_data, ...
  lmmcmaes_data, ...
  saacmes_data, ...
  dtscmaes_data};

datanames = {...
  'aDTS \epsilon_{min} Q_2 \epsilon_{max} Q_2', ... 'aDTS le50 he50 ur3', ...
  'aDTS \epsilon_{min} Q_2 \epsilon_{max} Q_3', ... 'aDTS le50 he75 ur3', ...
  'aDTS \epsilon_{min} Q_3 \epsilon_{max} Q_2', ... 'aDTS le75 he50 ur3', ...
  'aDTS \epsilon_{min} Q_3 \epsilon_{max} Q_3', ... 'aDTS le75 he75 ur3', ...
  'CMA-ES', ...
  'lmm-CMA-ES', ...
  '{}^{s*}ACMES', ...
  'DTS-CMA-ES'};

colors = [...
  adts_le50_he50_ur3_col; ...
  adts_le50_he75_ur3_col; ...
  adts_le75_he50_ur3_col; ...
  adts_le75_he75_ur3_col; ...
  cmaesCol; ...
  lmmCol; ...
  saacmesCol; ...
  dtsCol]/255;

testFunctions_2D  = [ 3,  4,  9, 12, 13, 14, ...
                     15, 16, 17, 21, 22, 23];
testFunctions_3D  = [ 3,  4, 11, 12, 13, 14, ...
                     16, 17, 18, 19, 21, 22];
testFunctions_5D  = [ 3,  8, 12, 13, 14, 15, ...
                     16, 17, 18, 19, 21, 22];
testFunctions_10D = [ 3,  4, 11, 12, 13, 14, ...
                     15, 17, 18, 19, 22, 24];
                  
if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% ADTS-CMA-ES settings comparison: lowErr, highErr, updateRate
% Ranking of function values in dimensions 2, 3, 5, 10.

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

set_datanames = {
  '$\lowErrMed$, $\highErrMed$, $\uRate = 0.3$',...
  '$\lowErrMed$, $\highErrMed$, $\uRate = 0.4$',...
  '$\lowErrMed$, $\highErrMed$, $\uRate = 0.5$',...
  '$\lowErrMed$, $\highErrQrt$, $\uRate = 0.3$',...
  '$\lowErrMed$, $\highErrQrt$, $\uRate = 0.4$',...
  '$\lowErrMed$, $\highErrQrt$, $\uRate = 0.5$',...
  '$\lowErrQrt$, $\highErrMed$, $\uRate = 0.3$',...
  '$\lowErrQrt$, $\highErrMed$, $\uRate = 0.4$',...
  '$\lowErrQrt$, $\highErrMed$, $\uRate = 0.5$',...
  '$\lowErrQrt$, $\highErrQrt$, $\uRate = 0.3$',...
  '$\lowErrQrt$, $\highErrQrt$, $\uRate = 0.4$',...
  '$\lowErrQrt$, $\highErrQrt$, $\uRate = 0.5$'};

% 2D

tableFunc = testFunctions_2D;
tableDims = 2;           
           
resultTable = fullfile(tableFolder, 'rankTable_2D.tex');
      
[table_2D, ranks] = rankingTable(set_data, 'DataNames', set_datanames, ...
                           'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                           'TableFuns', tableFunc, 'TableDims', tableDims,...
                           'Evaluations', [25, 50, 100, 200], ...
                           'ResultFile', resultTable, ...
                           'Mode', 'evaluations');
                         
% 3D

tableFunc = testFunctions_3D;
tableDims = 3;

resultTable = fullfile(tableFolder, 'rankTable_3D.tex');
      
[table_3D, ranks] = rankingTable(set_data, 'DataNames', set_datanames, ...
                           'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                           'TableFuns', tableFunc, 'TableDims', tableDims,...
                           'Evaluations', [25, 50, 100, 200], ...
                           'ResultFile', resultTable, ...
                           'Mode', 'evaluations');
                         
% 5D

tableFunc = testFunctions_5D;
tableDims = 5;

resultTable = fullfile(tableFolder, 'rankTable_5D.tex');
      
[table_5D, ranks] = rankingTable(set_data, 'DataNames', set_datanames, ...
                           'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                           'TableFuns', tableFunc, 'TableDims', tableDims,...
                           'Evaluations', [25, 50, 100, 200], ...
                           'ResultFile', resultTable, ...
                           'Mode', 'evaluations');
                         
% 10D

tableFunc = testFunctions_10D;
tableDims = 10;

resultTable = fullfile(tableFolder, 'rankTable_10D.tex');
      
[table_10D, ranks] = rankingTable(set_data, 'DataNames', set_datanames, ...
                           'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                           'TableFuns', tableFunc, 'TableDims', tableDims,...
                           'Evaluations', [25, 50, 100, 200], ...
                           'ResultFile', resultTable, ...
                           'Mode', 'evaluations');
                         
overall_table = table_2D(:, 1:4) + table_3D(:, 1:4) + table_5D(:, 1:4) + table_10D(:, 1:4);
                         
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
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% Algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Scaled function values of f1-f24 in dimension 10.

plotFuns = 1:24;
plotDims = 10;

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

%% Aggregated algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Aggregated  scaled function values in dimensions 5.

plotFuns = testFunctions_5D;
plotDims = [5, 10]; % both dimensions are generated because of legend split

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
                            
print2pdf(han{1}, pdfNames{1}, 1)

%% Aggregated algorithm comparison: CMA-ES, lmm-CMA-ES, saACMES, DTS-CMA-ES  
% Aggregated  scaled function values in dimensions 10.

plotFuns = testFunctions_10D;
plotDims = [5, 10]; % both dimensions are generated because of legend split

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
                            
print2pdf(han{2}, pdfNames{2}, 1)

%% Multiple comparison of algorithms with a statistical posthoc test.

close all

tableFunc = testFunctions_10D;
tableDims = 10;

resultDuelTable = fullfile(tableFolder, 'duelTable.tex');

datanames = {...
  '$\\lowErrMed$, $\\highErrMed$',...
  '$\\lowErrMed$, $\\highErrQrt$',...
  '$\\lowErrQrt$, $\\highErrMed$',...
  '$\\lowErrQrt$, $\\highErrQrt$',...
  'CMA-ES', 'lmm-CMA-ES', '\\saACMES', 'DTS-CMA-ES'};

[table, ranks] = duelTable(data, 'DataNames', datanames, ...
                            'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims, ...
                            'TableFuns', tableFunc, 'TableDims', tableDims, ...
                            'Evaluations', [1/3 1], ...
                            'ResultFile', resultDuelTable);

%% final clearing
close all
