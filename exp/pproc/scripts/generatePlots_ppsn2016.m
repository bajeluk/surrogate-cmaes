%% PPSN 2016 article plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values and graphs showing the speed up
% of GP DTS-CMA-ES.
% 
% Created for PPSN 2016 article.

%% load data

% path settings
exppath = fullfile('exp', 'experiments');

sd2_path = fullfile(exppath, 'exp_restrEC_04');
sd2_path20D = fullfile(exppath, 'exp_doubleEC_01_20D');
ei_poi_lcb_path = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb');
ei_poi_lcb_path20D = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb_20D');
cmaespath = fullfile(ei_poi_lcb_path, 'cmaes_results');
cmaespath20D = fullfile(ei_poi_lcb_path20D, 'cmaes_results');

saacmes_path = fullfile(exppath, 'saACMES');
smac_path = fullfile(exppath, 'SMAC');

% folder for results
plotResultsFolder = '/tmp';

% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10];

% loading data
[sd2_evals, sd2_settings] = dataReady(sd2_path, funcSet);
[ei_poi_lcb_evals, ei_poi_lcb_settings] = dataReady(ei_poi_lcb_path, funcSet);
cmaes_evals = dataReady(cmaespath, funcSet, 'cmaes');

funcSet.dims = 20;
[sd2_evals_20D, sd2_settings_20D] = dataReady(sd2_path20D, funcSet);
[ei_poi_lcb_evals_20D, ei_poi_lcb_settings_20D] = dataReady(ei_poi_lcb_path20D, funcSet);
cmaes_evals_20D = dataReady(cmaespath20D, funcSet, 'cmaes'); 

% funcSet.dims = [2, 3, 5, 10, 20];
% saacmes_evals = bbobDataReady(saacmes_path, funcSet);
% smac_evals = readSMACResults(smac_path, funcSet);

% finding data indexes
set.modelType = 'gp';
set.modelOpts.normalizeY = true;
set.evoControlRestrictedParam = 0.1;

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
sd2Id = getStructIndex(sd2_settings, set);
sd2Id20D = getStructIndex(sd2_settings_20D, set);

% concatenate data
eiData  = [ei_poi_lcb_evals(:, :, eiId),  ei_poi_lcb_evals_20D(:, :, eiId20D)];
poiData = [ei_poi_lcb_evals(:, :, poiId), ei_poi_lcb_evals_20D(:, :, poiId20D)];
lcbData = [ei_poi_lcb_evals(:, :, lcbId), ei_poi_lcb_evals_20D(:, :, lcbId20D)];
sd2Data = [sd2_evals(:, :, sd2Id), sd2_evals_20D(:, :, sd2Id20D)];
% saacmesData = saacmes_evals;
% smacData = smac_evals;
cmaesData = [cmaes_evals, cmaes_evals_20D];

%% f-values EI, PoI, lcb, sd2 comparison

close all

% color settings
cmaesCol = [22 22 138];
eiCol = [255 0 0];
poiCol = [255 215 0];
lcbCol = [208 32 144];
sd2Col = [0 0 0];
saacmesCol = [100 149 237];
smacCol = [116 172 66];

data = {eiData, ...
        poiData, ...
        lcbData, ...
        sd2Data, ...
        cmaesData};
%         smac_evals, ...
%         saacmes_evals, ...
%         cmaes_evals};

% data = {eiData, ...
%         ei_poi_lcb_evals(:, :, poiId), ...
%         ei_poi_lcb_evals(:, :, lcbId), ...
%         saacmes_evals, ...
%         cmaes_evals};
datanames = {'EI', 'poi', 'lcb', 'sd2', 'CMA-ES'};

% data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId),trans_evals(:,:,gp5TransId)};
% datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};

%   data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId)};
%   datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};


colors = [eiCol; poiCol; lcbCol; sd2Col; cmaesCol]/255;
% for i = 1:length(funcSet.BBfunc)
%   pdfNames{i} = fullfile(plotResultsFolder, ['f', num2str(funcSet.BBfunc(i))]);
% end

han = fValuesPlot(data, 'DataNames', datanames, 'DataDims', funcSet.dims, ...
                        'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                        'Statistic', 'median', 'AverageDims', false, ...
                        'Dependency', 'dim');
%   print2pdf(han,pdfNames,1)

%   drawGraph(gp_evals,cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [1,2,3,5,6,8];
%   h1 = drawComparison(trans_evals(:,:,1),trans_evals(:,:,2),cmaes_evals,'cmaes',funcSet);
%   
%   funcSet.BBfunc = [10,11,12,13,14,20,21];
%   h2 = drawComparison(gp_evals,rf_evals,cmaes_evals,'cmaes',funcSet);

%   pdfname = fullfile(plotResultsFolder,'speedUp');
% print2pdf(han, pdfNames, 1)

%% final clearing
close all