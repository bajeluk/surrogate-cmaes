%% GECC0 2021 article plots: Interaction between Model and Its Management in Surrogate-Assisted CMA Evolution Strategy
%
% *Paper abstract:*
% 
% Here we provide additional material which was not fully incorporated in
% the original paper.
%
% Created for GECCO 2021 article.

%% 

% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'gecco2021_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% needed function and dimension settings
funcSet.BBfunc = [1:24, 101:130, 207];
funcSet.dims = [2, 3, 4, 5, 10, 12, 16, 20];
maxEvals = 250;

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'gecco2021paper');
plotResultsFolder = fullfile(articleFolder, 'images');
tableFolder = fullfile(articleFolder, 'tex');
[~, ~] = mkdir(plotResultsFolder);
[~, ~] = mkdir(tableFolder);

% path settings
exppath = fullfile('exp', 'experiments', 'hanus', 'experiments');

% noiseless DT EC
nl_dts_gp_poi_path = fullfile(exppath, 'exp_doubleEC_26_1model');
nl_dts_gp_fval_path = fullfile(exppath, 'marek_doubleEC_gp_fvalues');
nl_dts_han_path = fullfile(exppath, 'marek_doubleEC_hansen');
nl_dts_lmm_path = fullfile(exppath, 'marek_doubleEC_lmm');
% noiseless linQuadEC
nl_lq_gp_poi_path = fullfile(exppath, 'marek_linQuadEC_gp_poi');
nl_lq_gp_fval_path = fullfile(exppath, 'marek_linQuadEC_gp');
nl_lq_han_path = fullfile(exppath, 'marek_linQuadEC_hansen');
nl_lq_lmm_path = fullfile(exppath, 'marek_linQuadEC_lmm');
% noiseless lmmEC
nl_lmm_gp_poi_path = fullfile(exppath, 'marek_lmmEC_gp_poi');
nl_lmm_gp_fval_path = fullfile(exppath, 'marek_lmmEC_gp');
nl_lmm_han_path = fullfile(exppath, 'marek_lmmEC_hansen');
nl_lmm_lmm_path = fullfile(exppath, 'marek_lmmEC_lmm');

% noisy DT EC
no_dts_gp_poi_path = fullfile(exppath, 'marek_noisy_doubleEC_gp');
no_dts_gp_fval_path = fullfile(exppath, 'marek_noisy_doubleEC_gp_fvalues');
no_dts_han_path = fullfile(exppath, 'marek_noisy_doubleEC_hansen');
no_dts_lmm_path = fullfile(exppath, 'marek_noisy_doubleEC_lmm');
% noisy linQuadEC
no_lq_gp_poi_path = fullfile(exppath, 'marek_noisy_linQuadEC_gp_poi');
no_lq_gp_fval_path = fullfile(exppath, 'marek_noisy_linQuadEC_gp');
no_lq_han_path = fullfile(exppath, 'marek_noisy_linQuadEC_hansen');
no_lq_lmm_path = fullfile(exppath, 'marek_noisy_linQuadEC_lmm');
% noisy lmmEC
no_lmm_gp_poi_path = fullfile(exppath, 'marek_noisy_lmmEC_gp_poi');
no_lmm_gp_fval_path = fullfile(exppath, 'marek_noisy_lmmEC_gp');
no_lmm_han_path = fullfile(exppath, 'marek_noisy_lmmEC_hansen');
no_lmm_lmm_path = fullfile(exppath, 'marek_noisy_lmmEC_lmm');

% newfunc207 DT EC
f207_dts_gp_poi_path = fullfile(exppath, 'marek_newfunc207_doubleEC_gp');
f207_dts_gp_fval_path = fullfile(exppath, 'marek_newfunc207_doubleEC_gp_fvalues');
f207_dts_han_path = fullfile(exppath, 'marek_newfunc207_doubleEC_hansen');
f207_dts_lmm_path = fullfile(exppath, 'marek_newfunc207_doubleEC_lmm');
% newfunc207 linQuadEC
f207_lq_gp_poi_path = fullfile(exppath, 'marek_newfunc207_linQuadEC_gp_poi');
f207_lq_gp_fval_path = fullfile(exppath, 'marek_newfunc207_linQuadEC_gp');
f207_lq_han_path = fullfile(exppath, 'marek_newfunc207_linQuadEC_hansen');
f207_lq_lmm_path = fullfile(exppath, 'marek_newfunc207_linQuadEC_lmm');
% newfunc207 lmmEC
f207_lmm_gp_poi_path = fullfile(exppath, 'marek_newfunc207_lmmEC_gp_poi');
f207_lmm_gp_fval_path = fullfile(exppath, 'marek_newfunc207_lmmEC_gp');
f207_lmm_han_path = fullfile(exppath, 'marek_newfunc207_lmmEC_hansen');
f207_lmm_lmm_path = fullfile(exppath, 'marek_newfunc207_lmmEC_lmm');

% load data
dataFolders = {nl_dts_gp_poi_path, ... noiseless DTS
               nl_dts_gp_fval_path, ...
               nl_dts_han_path, ...
               nl_dts_lmm_path, ...
               nl_lq_gp_poi_path, ... noiseless linQuadEC
               nl_lq_gp_fval_path, ...
               nl_lq_han_path, ...
               nl_lq_lmm_path, ...
               nl_lmm_gp_poi_path, ... noiseless lmm
               nl_lmm_gp_fval_path, ...
               nl_lmm_han_path, ...
               nl_lmm_lmm_path, ...
               no_dts_gp_poi_path, ... noisy DTS
               no_dts_gp_fval_path, ...
               no_dts_han_path, ...
               no_dts_lmm_path, ...
               no_lq_gp_poi_path, ... noisy linQuadEC
               no_lq_gp_fval_path, ...
               no_lq_han_path, ...
               no_lq_lmm_path, ...
               no_lmm_gp_poi_path, ... noisy lmmEC
               no_lmm_gp_fval_path, ...
               no_lmm_han_path, ...
               no_lmm_lmm_path, ...
               f207_dts_gp_poi_path, ... newfunc207 DT EC
               f207_dts_gp_fval_path, ...
               f207_dts_han_path, ...
               f207_dts_lmm_path, ...
               f207_lq_gp_poi_path, ... newfunc207 linQuadEC
               f207_lq_gp_fval_path, ...
               f207_lq_han_path, ...
               f207_lq_lmm_path, ...
               f207_lmm_gp_poi_path, ... newfunc207 lmmEC
               f207_lmm_gp_fval_path, ...
               f207_lmm_han_path, ...
               f207_lmm_lmm_path, ...
              };

[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings
clear findSet
% DTS EC
findSet.evoControl = 'doubletrained';
findSet.modelType = 'hansen';
dts_han_Id = getStructIndex(settings, findSet);
findSet.modelType = 'lmm';
dts_lmm_Id = getStructIndex(settings, findSet);
findSet.modelType = 'gp';
findSet.modelOpts.predictionType = 'poi';
dts_gp_poi_Id = getStructIndex(settings, findSet);
findSet.modelOpts.predictionType = 'fvalues';
dts_gp_fval_Id = getStructIndex(settings, findSet);

clear findSet
% linQuadEC
findSet.evoControl = 'linquad';
findSet.modelType = 'hansen';
lq_han_Id = getStructIndex(settings, findSet);
findSet.modelType = 'lmm';
lq_lmm_Id = getStructIndex(settings, findSet);
findSet.modelType = 'gp';
findSet.modelOpts.predictionType = 'poi';
lq_gp_poi_Id = getStructIndex(settings, findSet);
findSet.modelOpts.predictionType = 'fvalues';
lq_gp_fval_Id = getStructIndex(settings, findSet);

clear findSet
% lmm EC
findSet.evoControl = 'lmm';
findSet.modelType = 'hansen';
lmm_han_Id = getStructIndex(settings, findSet);
findSet.modelType = 'lmm';
lmm_lmm_Id = getStructIndex(settings, findSet);
findSet.modelType = 'gp';
findSet.modelOpts.predictionType = 'poi';
lmm_gp_poi_Id = getStructIndex(settings, findSet);
findSet.modelOpts.predictionType = 'fvalues';
lmm_gp_fval_Id = getStructIndex(settings, findSet);

% extract data
dts_gp_poi_data = evals(:, :, dts_gp_poi_Id(1));
dts_gp_poi_data_2 = evals(:, :, dts_gp_poi_Id(2));
% DTS GP PoI is composed of two sligtly different settings -> unify
use_Id_2 = cellfun(@(x, y) numel(x) < numel(y), evals(:, :, dts_gp_poi_Id(1)), evals(:, :, dts_gp_poi_Id(2)));
dts_gp_poi_data(use_Id_2) = dts_gp_poi_data_2(use_Id_2);
% the rest is simple
dts_gp_fval_data = evals(:, :, dts_gp_fval_Id);
dts_han_data = evals(:, :, dts_han_Id);
dts_lmm_data = evals(:, :, dts_lmm_Id);
lq_gp_poi_data = evals(:, :, lq_gp_poi_Id);
lq_gp_fval_data = evals(:, :, lq_gp_fval_Id);
lq_han_data = evals(:, :, lq_han_Id);
lq_lmm_data = evals(:, :, lq_lmm_Id);
lmm_gp_poi_data = evals(:, :, lmm_gp_poi_Id);
lmm_gp_fval_data = evals(:, :, lmm_gp_fval_Id);
lmm_han_data = evals(:, :, lmm_han_Id);
lmm_lmm_data = evals(:, :, lmm_lmm_Id);

% color settings
dts_gp_poi_Col  = [255, 165,   0];  % orange (#ffa500)
dts_gp_fval_Col = [255,   0,   0];  % light red (#ff0000)
dts_han_Col     = [255,   0, 255];  % magenta (#ff00ff)
dts_lmm_Col     = [154, 205,  50];  % yellow grass green (#9acd32)
lq_gp_poi_Col   = [  0,   0, 255];  % middle blue (#0000ff)
lq_gp_fval_Col  = [133, 55,  106];  % dark violet (#85376a)
lq_han_Col      = [  0,   0, 128];  % navy blue (#000080)
lq_lmm_Col      = [  0, 191, 191];  % light petroleum (#00bfbf)
lmm_gp_poi_Col  = [100, 149, 237];  % cornflower blue (#6495ed)
lmm_gp_fval_Col = [173, 255,  47];  % shining light green (#adff2f)
lmm_han_Col     = [ 12, 240, 248];  % light azure
lmm_lmm_Col     = [255, 192, 203];  % solomon pink (#ffc0cb)
% fminconCol   = getAlgColors(23); % 23=middle yellow
% fmincon_pureCol = [  0, 127,   0]; % dark forrest green (#007f00)

% marker settings
dts_gp_poi_Mark  = 'p';
dts_gp_fval_Mark = '>';
dts_han_Mark     = 'd';
dts_lmm_Mark     = '^';
lq_gp_poi_Mark   = 'v';
lq_gp_fval_Mark  = '<';
lq_han_Mark      = 'x';
lq_lmm_Mark      = '+';
lmm_gp_poi_Mark  = 'o';
lmm_gp_fval_Mark = 's';
lmm_han_Mark     = '<';
lmm_lmm_Mark     = 'd';

% load results of pairwise comparison of algorithms with Wilcoxon's test
% (done by Martin)
wcxResFile = fullfile('exp', 'experiments', 'gecco2021data', 'zbynek.mat');
wcx = load(wcxResFile);

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% EC and model combinations comparison graphs
% Function values across in dimension 10.

data = { ...
        dts_gp_poi_data, ...
        dts_gp_fval_data, ...
        dts_han_data, ...
        dts_lmm_data, ...
        lq_gp_poi_data, ...
        lq_gp_fval_data, ...
        lq_han_data, ...
        lq_lmm_data, ...
        lmm_gp_poi_data, ...
        lmm_gp_fval_data, ...
        lmm_han_data, ...
        lmm_lmm_data, ...
       };

datanames = { ...
    '$\mathrm{DT}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{PoI}$', ...
    '$\mathrm{DT}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{fval}$', ...
    '$\mathrm{DT}_\mathrm{EC}\!+\!\mathrm{lq}_\mathrm{M}$', ...
    '$\mathrm{DT}_\mathrm{EC}\!+\!\mathrm{lmm}_\mathrm{M}$', ...
    '$\mathrm{lq}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{PoI}$', ...
    '$\mathrm{lq}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{fval}$', ...
    '$\mathrm{lq}_\mathrm{EC}\!+\!\mathrm{lq}_\mathrm{M}$', ...
    '$\mathrm{lq}_\mathrm{EC}\!+\!\mathrm{lmm}_\mathrm{M}$', ...
    '$\mathrm{lmm}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{PoI}$', ...
    '$\mathrm{lmm}_\mathrm{EC}\!+\!\mathrm{GP}_\mathrm{M}^\mathrm{fval}$', ...
    '$\mathrm{lmm}_\mathrm{EC}\!+\mathrm{lq}_\mathrm{M}$', ...
    '$\mathrm{lmm}_\mathrm{EC}\!+\!\mathrm{lmm}_\mathrm{M}$' ...
    };

colors = [dts_gp_poi_Col; dts_gp_fval_Col; dts_han_Col; dts_lmm_Col; ...
          lq_gp_poi_Col;  lq_gp_fval_Col;  lq_han_Col;  lq_lmm_Col; ...
          lmm_gp_poi_Col; lmm_gp_fval_Col; lmm_han_Col; lmm_lmm_Col]/255;
markers = {dts_gp_poi_Mark; dts_gp_fval_Mark; dts_han_Mark; dts_lmm_Mark; ...
          lq_gp_poi_Mark;  lq_gp_fval_Mark;  lq_han_Mark;  lq_lmm_Mark; ...
          lmm_gp_poi_Mark; lmm_gp_fval_Mark; lmm_han_Mark; lmm_lmm_Mark};

plotDims = 10;
plotFuns = [2, 6, 9, 13, 17, 22, 107, 119, 207];

pdfNames = {};
for d = plotDims
  for f = plotFuns
    pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('alg_f%d_%dD', f, d));
  end
end

close all
han = relativeFValuesPlot(data, ...
                              'AggregateDims', false, ...
                              'Colors', colors, ...
                              'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, ...
                              'DataNames', datanames, ...
                              'FunctionNames', true, ...
                              'LegendOption', 'split', ...
                              'LineSpecification', '-', ...
                              'Markers', markers, ...
                              'MaxEval', maxEvals, ...
                              'MaxInstances', Inf, ...
                              'OneFigure', false, ...
                              'PlotDims', plotDims, ...
                              'PlotFuns', plotFuns(1:end-1), ...
                              'PlotGrid', [], ...
                              'Quantiles', true(1, 12), ...
                              'Statistic', 'quantile' ...
                       );

han(end+1) = relativeFValuesPlot(data, ...
                              'AggregateDims', false, ...
                              'Colors', colors, ...
                              'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, ...
                              'DataNames', datanames, ...
                              'FunctionNames', true, ...
                              'LegendOption', 'hide', ...
                              'LineSpecification', '-', ...
                              'Markers', markers, ...
                              'MaxEval', maxEvals, ...
                              'MaxInstances', Inf, ...
                              'OneFigure', false, ...
                              'PlotDims', plotDims, ...
                              'PlotFuns', 207, ...
                              'PlotGrid', [], ...
                              'Quantiles', true(1, 12), ...
                              'Statistic', 'quantile', ...
                              'TitleString', 'fsim Buoy placement simulation' ...
                       );

print2pdf(han, pdfNames, 1)

%close all

%% Multiple comparison of algorithms with a statistical posthoc test.
close all

% wcx

combnames = { ...
    '$\\ECM{\\dtEC}{\\gpfM}$', ...
    '$\\ECM{\\dtEC}{\\gppM}$', ...
    '$\\ECM{\\dtEC}{\\lqM}$', ...
    '$\\ECM{\\dtEC}{\\lmM}$', ...
    '$\\ECM{\\lqEC}{\\gpfM}$', ...
    '$\\ECM{\\lqEC}{\\gppM}$', ...
    '$\\ECM{\\lqEC}{\\lqM}$', ...
    '$\\ECM{\\lqEC}{\\lmM}$', ...
    '$\\ECM{\\lmEC}{\\gpfM}$', ...
    '$\\ECM{\\lmEC}{\\gppM}$', ...
    '$\\ECM{\\lmEC}{\\lqM}$', ...
    '$\\ECM{\\lmEC}{\\lmM}$' ...
    };
  
modelnames = {...
  '$\\gpfM$', ...
  '$\\gppM$', ...
  '$\\lqM$', ...
  '$\\lmM$' ...
  };

ecnames = {...
  '$\\dtEC$', ...
  '$\\lqEC$', ...
  '$\\lmEC$', ...
  };

% number of function evaluations ordering
evalBudgets = [250, 25, 50];%[1, 1/10, 1/5];
nFEorder = [3, 1];
allDims = [2, 3, 4, 5, 10, 12, 16, 20];

wcxFields = fieldnames(wcx);

% data loop
for wf = 1:numel(wcxFields)
  
  % actual data field
  dataField = wcxFields{wf};
  wcxCombFields = fieldnames(wcx.(dataField));
  
  % combination loop
  for cf = 1:numel(wcxCombFields)
    
    % actual combination field
    combField = wcxCombFields{cf};
    
    % set dimensions according to benchmarks
    switch dataField
      case 'OutputInformation'
        continue
      case 'Comparisons207'
        dataDims = [2, 4, 10, 12, 16, 20];
      case 'ComparisonsAll'
        continue
      otherwise
        dataDims = [2, 3, 5, 10, 20];
    end
    % set data captions according to benchmarks
    switch dataField
      % all noiseless
      case 'ComparisonsNoiseless'
        dataCaption = 'noiseless COCO benchmarks';
      % all noisy
      case 'ComparisonsNoisy'
        dataCaption = 'noisy COCO benchmarks';
      % f207
      case 'Comparisons207'
        dataCaption = 'simulation benchmark';
      otherwise
        continue
    end
    % set names and width according to data type
    if contains(combField, 'Combination')
      dataNames = combnames;
      oneColumn = false;
      printHeader = true;
      headColW = '1.9cm';
    elseif contains(combField, 'EvolutionControl')
      dataNames = ecnames;
      oneColumn = true;
      printHeader = false;
      headColW = '1.7cm';
    elseif contains(combField, 'Model')
      dataNames = modelnames;
      oneColumn = true;
      printHeader = false;
      headColW = '1.7cm';
    else
      continue
    end

    % combinations of dimensions
    if contains(combField, 'Summary')
      % final data settings
      shiftedArray = arrayfun(@(x) shiftdim(wcx.(dataField)(x).(combField), -1), nFEorder, 'Uni', false);
      actualData = cat(4, shiftedArray{:});
      % resulting file name
      resultTableName = sprintf('dt_%s_%s.tex', dataField, combField);
      resultDuelTable = fullfile(tableFolder, resultTableName);
      headerFileName = sprintf('dt_%s_%s_header.tex', dataField, combField);
      % print duel table
      wcxDuelTable(cat(3, actualData), ...
                                  'DataCaption', dataCaption, ...
                                  'DataDims', num2cell(dataDims), ...
                                  'DataNames', dataNames, ...
                                  'DefFile', fullfile(tableFolder, 'defFile.tex'), ...
                                  'Evaluations', evalBudgets(nFEorder), ...
                                  'HeadColW', headColW, ...
                                  'HeaderFile', fullfile(tableFolder, headerFileName), ...
                                  'OneColumn', oneColumn, ...
                                  'PrintHeader', printHeader, ...
                                  'ResultFile', resultDuelTable);
    end
    
  end
end

%% final clearing
close all
