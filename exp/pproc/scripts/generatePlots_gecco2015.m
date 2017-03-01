%% GECCO 2015 poster plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and graphs showing the speed up
% of GP and RF S-CMA-ES generation strategy.
% 
% Created for GECCO 2015 poster.

%% F-values

close all

% path settings
exppath = fullfile('exp', 'experiments');
gppath = fullfile(exppath, 'exp_geneEC_06_all', 'exp_geneEC_06');
% rfpath = {[gppath,'_rflite'], [gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
rfpath = {[gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
transPath10D = fullfile(exppath, 'exp_geneEC_08_10D');
transPath20D = fullfile(exppath, 'exp_geneEC_08_20D');
cmaespath = fullfile(gppath, 'cmaes_results');
plotResultsFolder = '/tmp';
% plotResultsFolder = fullfile('..', 'latex_scmaes', 'gecco2015poster', 'images');

% needed function and dimension settings
funcSet.BBfunc = [1, 2, 3, 5, 6, 8, 10, 11, 12, 13, 14, 20, 21];
funcSet.dims = [2, 5, 10];

% loading data
[trans_evals, trans_settings] = dataReady(transPath10D, funcSet);
[rf_evals, rf_settings] = dataReady(rfpath, funcSet);
[gp_evals, gp_settings] = dataReady(gppath, funcSet);
cmaes_evals = dataReady(cmaespath, funcSet);

% finding data indexes
set.modelType = 'rf';
set.evoControlModelGenerations = 1;
rf1TransId = getStructIndex(trans_settings,set);
rf1Id = getStructIndex(rf_settings,set);

set.evoControlModelGenerations = 5;
rf5TransId = getStructIndex(trans_settings,set);

set.modelType = 'gp';
set.evoControlSampleRange = 1;
gp5TransId = getStructIndex(trans_settings,set);

set.evoControlModelGenerations = 1;
gp1TransId = getStructIndex(trans_settings,set);
gp1Id = getStructIndex(gp_settings,set);
set.evoControlModelGenerations = 3;
gp3Id = getStructIndex(gp_settings,set);

% color settings
CMAESCol = [22 22 138];
GP1TransCol = [255 0 0];
GP5TransCol = [255 215 0];
RF1TransCol = [208 32 144];
RF5TransCol = [0 0 0];
GP3Col = [100 149 237];
RF1Col = [116 172 66];

%   data = {trans_evals(:,:,rf1TransId),trans_evals(:,:,gp1TransId),trans_evals(:,:,rf5TransId),trans_evals(:,:,gp5TransId),cmaes_evals};
%   datanames = {'RF1T','GP1T','RF5T','GP5T','CMA-ES'};

data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId)};
datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};

%   data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId)};
%   datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};


colors = [CMAESCol; GP3Col; GP5TransCol; RF1Col; RF1TransCol]/255;
for i = 1:length(funcSet.BBfunc)
  pdfNames{i} = fullfile(plotResultsFolder, ['f', num2str(funcSet.BBfunc(i))]);
end

han = fValuesPlot(data, 'DataNames', datanames, 'DataDims', funcSet.dims, ...
                        'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                        'AverageDims', true);
%   print2pdf(han,pdfNames,1)

%   drawGraph(gp_evals,cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [1,2,3,5,6,8];
%   h1 = drawComparison(trans_evals(:,:,1),trans_evals(:,:,2),cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [10,11,12,13,14,20,21];
%   h2 = drawComparison(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);

%   pdfname = fullfile(plotResultsFolder,'speedUp');
print2pdf(han, pdfNames, 1)

%% Speed Up
% The speed up is considered between the CMA-ES and S-CMA-ES generation 
% strategy using GP, or RF.

% close all

exppath = fullfile('exp', 'experiments');
gppath = fullfile(exppath, 'exp_geneEC_06_all', 'exp_geneEC_06');
% rfpath = {[gppath,'_rflite'], [gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
rfpath = {[gppath,'_rf5_2'], [gppath,'_rf5_3'], [gppath,'_rf5_4']};
cmaespath = fullfile(gppath, 'cmaes_results');
plotResultsFolder = '/tmp';

funcSet.BBfunc = [1,2,3,5,6,8,10,11,12,13,14,20,21];
funcSet.BBfuncInv = inverseIndex(funcSet.BBfunc);
funcSet.dims = [2,5,10];
funcSet.dimsInv = inverseIndex(funcSet.dims);

rf_evals = dataReady(rfpath, funcSet, 1, 'rflite');
gp_evals = dataReady(gppath, funcSet, 6);
cmaes_evals = dataReady(cmaespath, funcSet, 1, 'cmaes');

hRFall = speedUpPlot(rf_evals, cmaes_evals, 'cmaes', funcSet);
pdfname = fullfile(plotResultsFolder, 'speedUpRF');
print2pdf(hRFall, pdfname, 1)

%   speedUpPlot(gp_evals,cmaes_evals,'cmaes',funcSet);

funcSet.BBfunc = [1,2,3,5,6,8];
h1 = speedUpPlotCompare(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);

funcSet.BBfunc = [10,11,12,14,20,21];
h2 = speedUpPlotCompare(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);

pdfname = fullfile(plotResultsFolder, 'speedUp');
print2pdf([h1, h2], {[pdfname, 'A.pdf'], [pdfname, 'B.pdf']}, 1)