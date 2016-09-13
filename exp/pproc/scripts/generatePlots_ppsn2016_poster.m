%% PPSN 2016 poster plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and graphs showing the speed up
% of GP DTS-CMA-ES.
% 
% Created for PPSN 2016 poster.

%% load data

% checkout file containing all loaded data
if ispc
  osTmp = fullfile('exp', 'pproc', 'scripts', 'tmp');
  if ~exist(osTmp, 'dir')
    mkdir(osTmp)
  end
else
  osTmp = '/tmp';
end
tmpFName = fullfile(osTmp, 'ppsndata.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'ppsn2016poster');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');

% path settings
exppath = fullfile('exp', 'experiments');

sd2_r10_20_path = fullfile(exppath, 'exp_restrEC_04');
sd2_r05_40_path = fullfile(exppath, 'exp_doubleEC_01_restr05_40');
sd2_r05_2pop_path = fullfile(exppath, 'exp_restrEC_04_2pop');
sd2_r10_2pop_path = fullfile(exppath, 'exp_doubleEC_01_2pop');
sd2_r20_40_2pop_path = fullfile(exppath, 'exp_doubleEC_01_2pop_restr20_40');
sd2_path20D = fullfile(exppath, 'exp_doubleEC_01_20D');

ei_poi_lcb_path = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb');
ei_poi_lcb_path20D = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb_20D');

gen_path = fullfile(exppath, 'exp_geneEC_10');
gen_path20D = fullfile(exppath, 'exp_geneEC_10_20D');

cmaespath = fullfile(gen_path, 'cmaes_results');
cmaespath20D = fullfile(ei_poi_lcb_path20D, 'cmaes_results');

saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
smac_path = fullfile(exppath, 'SMAC');

% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10];

% loading data
[sd2_r10_20_evals, sd2_r10_20_settings] = dataReady(sd2_r10_20_path, funcSet);
[sd2_r05_40_evals, sd2_r05_40_settings] = dataReady(sd2_r05_40_path, funcSet);
[sd2_r05_2pop_evals, sd2_r05_2pop_settings] = dataReady(sd2_r05_2pop_path, funcSet);
[sd2_r10_2pop_evals, sd2_r10_2pop_settings] = dataReady(sd2_r10_2pop_path, funcSet);
[sd2_r20_40_2pop_evals, sd2_r20_40_2pop_settings] = dataReady(sd2_r20_40_2pop_path, funcSet);

[ei_poi_lcb_evals, ei_poi_lcb_settings] = dataReady(ei_poi_lcb_path, funcSet);

[gen_evals, gen_settings] = dataReady(gen_path, funcSet);

% concatenate cmaes 2-10D
cmaes_evals = dataReady(cmaespath, funcSet);
for f = 1:size(cmaes_evals, 1)
  for d = 1:size(cmaes_evals, 2)
    notEmptyId = ~cellfun(@isempty, cmaes_evals(f, d, :));
    cmaes_evals{f,d,1} = cmaes_evals{f, d, notEmptyId};
  end
end
cmaes_evals = cmaes_evals(:, :, 1);

funcSet.dims = 20;
[sd2_evals_20D, sd2_settings_20D] = dataReady(sd2_path20D, funcSet);
[ei_poi_lcb_evals_20D, ei_poi_lcb_settings_20D] = dataReady(ei_poi_lcb_path20D, funcSet);
[gen_evals_20D, gen_settings_20D] = dataReady(gen_path20D, funcSet);
% Uncomment if necessary:
% This is a hack due to distributed and merged part of 20D experiment:
% if (length(gen_settings_20D) > 4)
%   gen_settings_20D(1:4) = gen_settings_20D((end-3):end);
%   gen_settings_20D(5:end) = [];
% end

% concatenate cmaes 20D
cmaes_evals_20D = dataReady(cmaespath20D, funcSet);
for f = 1:size(cmaes_evals_20D, 1)
  for d = 1:size(cmaes_evals_20D, 2)
    notEmptyId = ~cellfun(@isempty, cmaes_evals_20D(f, d, :));
    cmaes_evals_20D{f,d,1} = cmaes_evals_20D{f, d, notEmptyId};
  end
end
cmaes_evals_20D = cmaes_evals_20D(:, :, 1);

funcSet.dims = [2, 3, 5, 10, 20];
saacmes_evals = bbobDataReady(saacmes_path, funcSet);
smac_evals = bbobDataReady(smac_path, funcSet);

% finding data indexes
clear set
set.modelType = 'gp';
set.modelOpts.normalizeY = true;
set.evoControlModelGenerations = 5;
genId = getStructIndex(gen_settings, set);
genId20D = getStructIndex(gen_settings_20D, set);

set = rmfield(set, 'evoControlModelGenerations');
set.evoControlRestrictedParam = 0.1;
set.PopSize = '(4 + floor(3*log(N)))';

set.modelOpts.predictionType = 'ei';
eiId = getStructIndex(ei_poi_lcb_settings, set);
eiId20D = getStructIndex(ei_poi_lcb_settings_20D, set);

set.modelOpts.predictionType = 'poi';
poiId = getStructIndex(ei_poi_lcb_settings, set);
poiId20D = getStructIndex(ei_poi_lcb_settings_20D, set);

set.modelOpts.predictionType = 'lcb';
lcbId = getStructIndex(ei_poi_lcb_settings, set);
lcbId20D = getStructIndex(ei_poi_lcb_settings_20D, set);

set.modelOpts.predictionType = 'sd2';
set.evoControlRestrictedParam = 0.05;
sd2_r05_Id = getStructIndex(sd2_r05_40_settings, set);
sd2_r05_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.1;
sd2_r10_Id = getStructIndex(sd2_r10_20_settings, set);
sd2_r10_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.2;
sd2_r20_Id = getStructIndex(sd2_r10_20_settings, set);
sd2_r20_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.4;
sd2_r40_Id = getStructIndex(sd2_r05_40_settings, set);
sd2_r40_Id20D = getStructIndex(sd2_settings_20D, set);

set.PopSize = '(8 + floor(6*log(N)))';
set.evoControlRestrictedParam = 0.05;
sd2_r05_2pop_Id = getStructIndex(sd2_r05_2pop_settings, set);
sd2_r05_2pop_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.1;
sd2_r10_2pop_Id = getStructIndex(sd2_r10_2pop_settings, set);
sd2_r10_2pop_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.2;
sd2_r20_2pop_Id = getStructIndex(sd2_r20_40_2pop_settings, set);
sd2_r20_2pop_Id20D = getStructIndex(sd2_settings_20D, set);

set.evoControlRestrictedParam = 0.4;
sd2_r40_2pop_Id = getStructIndex(sd2_r20_40_2pop_settings, set);
sd2_r40_2pop_Id20D = getStructIndex(sd2_settings_20D, set);

% concatenate data
eiData  = [ei_poi_lcb_evals(:, :, eiId),  ei_poi_lcb_evals_20D(:, :, eiId20D)];
poiData = [ei_poi_lcb_evals(:, :, poiId), ei_poi_lcb_evals_20D(:, :, poiId20D)];
lcbData = [ei_poi_lcb_evals(:, :, lcbId), ei_poi_lcb_evals_20D(:, :, lcbId20D)];

sd2Data_05 = [sd2_r05_40_evals(:, :, sd2_r05_Id), sd2_evals_20D(:, :, sd2_r05_Id20D)];
sd2Data_10 = [sd2_r10_20_evals(:, :, sd2_r10_Id), sd2_evals_20D(:, :, sd2_r10_Id20D)];
sd2Data_20 = [sd2_r10_20_evals(:, :, sd2_r20_Id), sd2_evals_20D(:, :, sd2_r20_Id20D)];
sd2Data_40 = [sd2_r05_40_evals(:, :, sd2_r40_Id), sd2_evals_20D(:, :, sd2_r40_Id20D)];

sd2Data_05_2pop = [sd2_r05_2pop_evals(:, :, sd2_r05_2pop_Id), sd2_evals_20D(:, :, sd2_r05_2pop_Id20D)];
sd2Data_10_2pop = [sd2_r10_2pop_evals(:, :, sd2_r10_2pop_Id), sd2_evals_20D(:, :, sd2_r10_2pop_Id20D)];
sd2Data_20_2pop = [sd2_r20_40_2pop_evals(:, :, sd2_r20_2pop_Id), sd2_evals_20D(:, :, sd2_r20_2pop_Id20D)];
sd2Data_40_2pop = [sd2_r20_40_2pop_evals(:, :, sd2_r40_2pop_Id), sd2_evals_20D(:, :, sd2_r40_2pop_Id20D)];

saacmesData = saacmes_evals;
smacData = smac_evals;
cmaesData = [cmaes_evals(:, :, 1) , cmaes_evals_20D(:, :, 1)];
genData = [gen_evals(:, :, genId), gen_evals_20D(:, :, genId20D)];

% color settings
cmaesCol = [22 22 138];

eiCol = [255 0 0];
poiCol = [255 215 0];
lcbCol = [208 32 144];
sd2Col = [0 0 0];
saacmesCol = [100 149 237];
smacCol = [255, 155, 0];
genCol = [178,34,34];

sd2Col_05 = [200 170 39];
sd2Col_10 = sd2Col;
sd2Col_20 = [148 0 211];
sd2Col_40 = [255 20 147];
sd2Col_05_2pop = [154 205 50];
sd2Col_10_2pop = [34,139,34];
sd2Col_20_2pop = [0,128,128];
sd2Col_40_2pop = [70,130,180];

% evaluation target settings
defTargets = floor(power(20, linspace(1, log(250)/log(20), 25)));

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Used Output

%% Criterion comparison: EI, PoI, lcb, sd2
% Aggregation of function values across dimensions  5, 20.

data = {eiData, ...
        poiData, ...
        lcbData, ...
        sd2Data_10, ...
        cmaesData};

datanames = {'EI', 'PoI', 'LCB', 's^2', 'CMA-ES'};

colors = [eiCol; poiCol; lcbCol; sd2Col; cmaesCol]/255;

plotDims = [10, 20];

clear pdfNames
pdfNames = arrayfun(@(x) fullfile(plotResultsFolder, ['crit', num2str(x), 'D']), plotDims, 'UniformOutput', false);

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LegendOption', 'show', 'MaxEval', 100);
                            
print2pdf(han, pdfNames, 1)

%% Population size comparison: default, 2*default
% Aggregation of function values across dimensions 10, 20 with splitted
% legend

data = {sd2Data_05, ...
        sd2Data_10, ...
        sd2Data_20, ...
        sd2Data_40, ...
        sd2Data_05_2pop, ...
        sd2Data_10_2pop, ...
        sd2Data_20_2pop, ...
        sd2Data_40_2pop, ...
        cmaesData};

datanames = {'0.05 1pop', '0.1  1pop', '0.2  1pop', '0.4  1pop', ...
             '0.05 2pop', '0.1  2pop', '0.2  2pop', '0.4  2pop', ...
             'CMA-ES'};

colors = [sd2Col_05; sd2Col_10; sd2Col_20; sd2Col_40; ...
          sd2Col_05_2pop; sd2Col_10_2pop; sd2Col_20_2pop; sd2Col_40_2pop; ...
          cmaesCol]/255;
        
plotDims = [10, 20];

clear pdfNames
pdfNames = arrayfun(@(x) fullfile(plotResultsFolder, ['pop', num2str(x), 'D']), plotDims, 'UniformOutput', false);

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', true, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LegendOption', 'split', 'MaxEval', 100);
                            
% print2pdf(han, pdfNames, 1)

%% Algorithm comparison: DTS-CMA-ES, S-CMA-ES, saACMES, SMAC, CMA-ES
% Aggregation of function values across dimensions 5, 10, 20.

data = {sd2Data_10, ...
        sd2Data_05_2pop, ...
        genData, ...
        saacmesData, ...
        smacData, ...
        cmaesData};

datanames = {'DTS 0.1  1pop', 'DTS 0.05 2pop', 'S-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'SMAC', 'CMA-ES'};

colors = [sd2Col_10; sd2Col_05_2pop; genCol; saacmesCol; smacCol; cmaesCol]/255;

plotDims = [2, 5, 10, 20];

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', true, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpec', {'-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'split', 'MaxEval', 100,...
                              'LineWidth', [2, 2, 1, 1, 1, 1]);

print2pdf(han, fullfile(plotResultsFolder, 'alg2_5_10_20D'), 1)

plotDims = [10, 20];

clear pdfNames
pdfNames = {fullfile(plotResultsFolder, 'alg10D'), ...
            fullfile(plotResultsFolder, 'alg20D')};

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', funcSet.BBfunc, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', true, ...
                              'LineSpec', {'-', '-', '-', '-', '-', '-'}, ...
                              'LegendOption', 'hide', 'MaxEval', 100,...
                              'LineWidth', [2, 2, 1, 1, 1, 1]);
                            
print2pdf(han, pdfNames, 1)

%% Algorithm comparison: DTS-CMA-ES, S-CMA-ES, saACMES, SMAC, CMA-ES
% Aggregation of dimension values across all functions.

data = {sd2Data_10, ...
        sd2Data_05_2pop, ...
        genData, ...
        saacmesData, ...
        smacData, ...
        cmaesData};

datanames = {'DTS 0.1  1pop', 'DTS 0.05 2pop', 'S-CMA-ES', 'BIPOP-{}^{s*}ACMES-k', 'SMAC', 'CMA-ES'};

colors = [sd2Col_10; sd2Col_05_2pop; genCol; saacmesCol; smacCol; cmaesCol]/255;

plotFuns = 1:24;
plotDims = [10];

clear pdfNames
pdfNames = arrayfun(@(x) fullfile(plotResultsFolder, ['alg_f', num2str(x)]), plotFuns, 'UniformOutput', false);

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, 'AggregateFuns', false, ...
                              'LineSpec', {'-', '-', '-', '-', '-', '-'}, ...
                              'LineWidth', [2, 2, 1, 1, 1, 1], ...
                              'LegendOption', 'split', 'MaxEval', 100, ...
                              'FunctionNames', true);

                              
                            
print2pdf(han, pdfNames, 1)

%% final clearing
close all
