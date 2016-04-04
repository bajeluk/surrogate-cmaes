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
% sd2_path20D = fullfile(exppath, 'exp_doubleEC_01_20D');
ei_poi_lcb_path = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb');
% ei_poi_lcb_path20D = fullfile(exppath, 'exp_doubleEC_01_ei_poi_lcb_20D');
cmaespath = fullfile(ei_poi_lcb_path, 'cmaes_results');

plotResultsFolder = '/tmp';

% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10];

% loading data
[sd2_evals, sd2_settings] = dataReady(sd2_path, funcSet, 4);
[ei_poi_lcb_evals, ei_poi_lcb_settings] = dataReady(ei_poi_lcb_path, funcSet, 3);
% [ei_poi_lcb20D_evals, ei_poi_lcb20D_settings] = dataReady(ei_poi_lcb_path20D, funcSet, 3);
cmaes_evals = dataReady(cmaespath, funcSet, 1, 'cmaes');

%% f-values EI, PoI, lcb, sd2 comparison

close all

% finding data indexes
set.modelType = 'gp';
set.modelOpts.normalizeY = true;
set.evoControlRestrictedParam = 0.1;

set.modelOpts.predictionType = 'ei';
eiId = getStructIndex(ei_poi_lcb_settings, set);

set.modelOpts.predictionType = 'poi';
poiId = getStructIndex(ei_poi_lcb_settings, set);

set.modelOpts.predictionType = 'lcb';
lcbId = getStructIndex(ei_poi_lcb_settings, set);

set.modelOpts.predictionType = 'sd2';
sd2Id = getStructIndex(sd2_settings, set);

% set.evoControlModelGenerations = 5;
% rf5TransId = getStructIndex(trans_settings,set);
% 
% set.modelType = 'gp';
% set.evoControlSampleRange = 1;
% gp5TransId = getStructIndex(trans_settings,set);
% 
% set.evoControlModelGenerations = 1;
% gp1TransId = getStructIndex(trans_settings,set);
% gp1Id = getStructIndex(gp_settings,set);
% set.evoControlModelGenerations = 3;
% gp3Id = getStructIndex(gp_settings,set);
% 

% color settings
CMAESCol = [22 22 138];
eiCol = [255 0 0];
poiCol = [255 215 0];
lcbCol = [208 32 144];
sd2Col = [0 0 0];
% GP3Col = [100 149 237];
% RF1Col = [116 172 66];

data = {ei_poi_lcb_evals(:, :, eiId), ...
        ei_poi_lcb_evals(:, :, poiId), ...
        ei_poi_lcb_evals(:, :, lcbId), ...
        sd2_evals(:, :, sd2Id), ...
        cmaes_evals};
datanames = {'EI', 'PoI', 'lcb', 'sd2', 'CMA-ES'};

% data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId),trans_evals(:,:,gp5TransId)};
% datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};

%   data = {cmaes_evals,gp_evals(:,:,gp3Id),trans_evals(:,:,gp5TransId),rf_evals(:,:,rf1Id),trans_evals(:,:,rf1TransId)};
%   datanames = {'CMA-ES','GP3','GP5-trans','RF1','RF1-trans'};


colors = [CMAESCol; eiCol; poiCol; lcbCol; sd2Col]/255;
% for i = 1:length(funcSet.BBfunc)
%   pdfNames{i} = fullfile(plotResultsFolder, ['f', num2str(funcSet.BBfunc(i))]);
% end

han = fValuesPlot(data, datanames, funcSet, funcSet.dims, funcSet.BBfunc, colors);
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

