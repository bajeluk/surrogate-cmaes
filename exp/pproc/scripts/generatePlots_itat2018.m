%% ITAT 2018 plots
% Script for making graphs showing the dependence of minimal function
% values on the number of function values of compared algorithms.
% 
% Created for ITAT 2018 workshop article.

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'itat2018_data.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else
  
% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10];
maxEvals = 250;
  
% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('surrogate-cmaes')), 'latex_scmaes', 'itat2018paper');
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

% gen_path = fullfile(exppath, 'exp_geneEC_10');
% gen_path20D = fullfile(exppath, 'exp_geneEC_10_20D');
dts_path = fullfile(exppath, 'DTS-CMA-ES_05_2pop');
dts_rf_path = fullfile(exppath, 'exp_doubleEC_rf_s01');
% maes_path = fullfile(exppath, 'exp_maesEC_14_2_10_cmaes_20D');

cmaes_path = fullfile(exppath, 'CMA-ES');
% saacmes_path = fullfile(exppath, 'BIPOP-saACM-k');
lmm_path = fullfile(exppath, 'lmm-CMA-ES');

% load data
dataFolders = {dts_path; ...
               dts_rf_path; ...
               cmaes_path; ...
               lmm_path};
             
[evals, settings] = catEvalSet(dataFolders, funcSet);

% find ids in settings
clear findSet
findSet.modelType = 'forest';
% restrictedParam = 0.05
findSet.evoControlRestrictedParam = 0.05;
findSet.modelOpts.tree_splitFunc = @AxisSplit;
rf_05_axisId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @GaussianSplit;
rf_05_gaussId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @HillClimbingObliqueSplit;
rf_05_hcId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @PairObliqueSplit;
rf_05_pairId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @ResidualObliqueSplit ;
rf_05_resId = getStructIndex(settings, findSet);

% restrictedParam = 0.10
findSet.evoControlRestrictedParam = 0.10;
findSet.modelOpts.tree_splitFunc = @AxisSplit;
rf_10_axisId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @GaussianSplit;
rf_10_gaussId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @HillClimbingObliqueSplit;
rf_10_hcId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @PairObliqueSplit;
rf_10_pairId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @ResidualObliqueSplit ;
rf_10_resId = getStructIndex(settings, findSet);

% restrictedParam = 0.20
findSet.evoControlRestrictedParam = 0.20;
findSet.modelOpts.tree_splitFunc = @AxisSplit;
rf_20_axisId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @GaussianSplit;
rf_20_gaussId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @HillClimbingObliqueSplit;
rf_20_hcId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @PairObliqueSplit;
rf_20_pairId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @ResidualObliqueSplit ;
rf_20_resId = getStructIndex(settings, findSet);

% restrictedParam = 0.40
findSet.evoControlRestrictedParam = 0.40;
findSet.modelOpts.tree_splitFunc = @AxisSplit;
rf_40_axisId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @GaussianSplit;
rf_40_gaussId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @HillClimbingObliqueSplit;
rf_40_hcId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @PairObliqueSplit;
rf_40_pairId = getStructIndex(settings, findSet);
findSet.modelOpts.tree_splitFunc = @ResidualObliqueSplit ;
rf_40_resId = getStructIndex(settings, findSet);

% reference algorithms Ids
clear findSet
findSet.algName = 'CMA-ES';
cma_Id = getStructIndex(settings, findSet);
% findSet.algName = 'BIPOP-saACM-k';
% saacm_Id = getStructIndex(settings, findSet);
findSet.algName = 'lmm-CMA-ES';
lmm_Id = getStructIndex(settings, findSet);
findSet.algName = 'DTS-CMA-ES_05_2pop';
dts_Id = getStructIndex(settings, findSet);
             
% extract data
rf_05_axis_data = evals(:, :, rf_05_axisId);
rf_05_gauss_data = evals(:, :, rf_05_gaussId);
rf_05_hc_data = evals(:, :, rf_05_hcId);
rf_05_pair_data = evals(:, :, rf_05_pairId);
rf_05_res_data = evals(:, :, rf_05_resId);

rf_10_axis_data = evals(:, :, rf_10_axisId);
rf_10_gauss_data = evals(:, :, rf_10_gaussId);
rf_10_hc_data = evals(:, :, rf_10_hcId);
rf_10_pair_data = evals(:, :, rf_10_pairId);
rf_10_res_data = evals(:, :, rf_10_resId);

rf_20_axis_data = evals(:, :, rf_20_axisId);
rf_20_gauss_data = evals(:, :, rf_20_gaussId);
rf_20_hc_data = evals(:, :, rf_20_hcId);
rf_20_pair_data = evals(:, :, rf_20_pairId);
rf_20_res_data = evals(:, :, rf_20_resId);

rf_40_axis_data = evals(:, :, rf_40_axisId);
rf_40_gauss_data = evals(:, :, rf_40_gaussId);
rf_40_hc_data = evals(:, :, rf_40_hcId);
rf_40_pair_data = evals(:, :, rf_40_pairId);
rf_40_res_data = evals(:, :, rf_40_resId);

cmaes_data     = evals(:, :, cma_Id);
% saacmes_data   = evals(:, :, saacm_Id);
dtscmaes_data  = evals(:, :, dts_Id);
lmmcmaes_data  = evals(:, :, lmm_Id);

% color settings
rf_05_axisCol = getAlgColors(1);
rf_05_gaussCol = getAlgColors(2);
rf_05_hcCol = getAlgColors(3);
rf_05_pairCol = getAlgColors(4);
rf_05_resCol = getAlgColors(5);

rf_10_axisCol = getAlgColors(6);
rf_10_gaussCol = getAlgColors(7);
rf_10_hcCol = getAlgColors(8);
rf_10_pairCol = getAlgColors(9);
rf_10_resCol = getAlgColors(10);

rf_20_axisCol = getAlgColors(11);
rf_20_gaussCol = getAlgColors(12);
rf_20_hcCol = getAlgColors(13);
rf_20_pairCol = getAlgColors(14);
rf_20_resCol = getAlgColors(15);

rf_40_axisCol = getAlgColors(16);
rf_40_gaussCol = getAlgColors(17);
rf_40_hcCol = getAlgColors(18);
rf_40_pairCol = getAlgColors(19);
rf_40_resCol = getAlgColors(20);

cmaesCol     = getAlgColors('cmaes');
% saacmesCol   = getAlgColors('saacmes');
dtsCol       = getAlgColors('dtscmaes');
lmmCol       = getAlgColors('lmmcmaes');

axisCol   = [255, 165,   0];  % orange (#ffa500)
gaussCol  = [255,   0,   0];  % light red (#ff0000)
hcCol     = [255,   0, 255];  % magenta (#ff00ff)
pairCol   = [  0,   0, 255];  % middle blue (#0000ff)
resCol    = [133, 55,  106];  % dark violet (#85376a)

% marker settings
axisMark  = 'p';
gaussMark = '>';
hcMark    = 'd';
pairMark  = '^';
resMark   = 'v';

rf05Mark  = '<';
rf10Mark  = 'd';
rf20Mark  = '>';
rf40Mark  = 'p';


% sapeoMark       = '<';
cmaesMark = 'x';
% cmaes2popMark   = '+';
dtsMark   = 'o';
lmmMark   = 's';

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% Split function comparison: Axis, Gauss, HC, Pair, Residual, CMA-ES
% Aggregation of function values across dimensions 2, 3, 5, 10.

axis_data = cellfun(@(x1, x2, x3, x4) [x1, x2, x3, x4], ...
  rf_05_axis_data, rf_10_axis_data, rf_20_axis_data, rf_40_axis_data, ...
  'UniformOutput', false);
gauss_data = cellfun(@(x1, x2, x3, x4) [x1, x2, x3, x4], ...
  rf_05_gauss_data, rf_10_gauss_data, rf_20_gauss_data, rf_40_gauss_data, ...
  'UniformOutput', false);
hc_data = cellfun(@(x1, x2, x3, x4) [x1, x2, x3, x4], ...
  rf_05_hc_data, rf_10_hc_data, rf_20_hc_data, rf_40_hc_data, ...
  'UniformOutput', false);
pair_data = cellfun(@(x1, x2, x3, x4) [x1, x2, x3, x4], ...
  rf_05_pair_data, rf_10_pair_data, rf_20_pair_data, rf_40_pair_data, ...
  'UniformOutput', false);
res_data = cellfun(@(x1, x2, x3, x4) [x1, x2, x3, x4], ...
  rf_05_res_data, rf_10_res_data, rf_20_res_data, rf_40_res_data, ...
  'UniformOutput', false);

data = { ...
        axis_data, ...
        gauss_data, ...
        hc_data, ...
        pair_data, ...
        res_data, ...
        cmaes_data, ...
        dtscmaes_data
        };
      
datanames = { ...
    'CART', ...
    'SECRET', ...
    'OC1', ...
    'PAIR', ...
    'SUPPORT', ...
    'CMA-ES', ...
    'DTS-CMA-ES' ...
    % 'adaptive DTS-CMA-ES' ...
    };
  
colors = [axisCol; gaussCol; hcCol; pairCol; resCol; cmaesCol; dtsCol]/255; 
markers = {axisMark, gaussMark, hcMark, pairMark, resMark, cmaesMark, dtsMark};

plotFuns = [1:24];
plotDims = [2, 3, 5];

clear pdfNames
pdfNames = {};
for d = plotDims
  pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('split_%dD', d));
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, ...
                              'Quantiles', false(1, 7), ...
                              'AggregateFuns', true, ...
                              'LineSpecification', '-', ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'Markers', markers, ...
                              'PlotGrid', [2, 2], ...
                              'MaxInstances', Inf, ...
                              'FunctionNames', true);
                            
print2pdf(han, pdfNames, 1)

%% Number of originally-evaluated points comparison
% Aggregation of function values across dimensions 2, 3, 5, 10.

rf_05_data = cellfun(@(x1, x2, x3, x4, x5) [x1, x2, x3, x4, x5], ...
  rf_05_axis_data, rf_05_gauss_data, rf_05_hc_data, rf_05_pair_data, rf_05_res_data, ...
  'UniformOutput', false);
rf_10_data = cellfun(@(x1, x2, x3, x4, x5) [x1, x2, x3, x4, x5], ...
  rf_10_axis_data, rf_10_gauss_data, rf_10_hc_data, rf_10_pair_data, rf_10_res_data, ...
  'UniformOutput', false);
rf_20_data = cellfun(@(x1, x2, x3, x4, x5) [x1, x2, x3, x4, x5], ...
  rf_20_axis_data, rf_20_gauss_data, rf_20_hc_data, rf_20_pair_data, rf_20_res_data, ...
  'UniformOutput', false);
rf_40_data = cellfun(@(x1, x2, x3, x4, x5) [x1, x2, x3, x4, x5], ...
  rf_40_axis_data, rf_40_gauss_data, rf_40_hc_data, rf_40_pair_data, rf_40_res_data, ...
  'UniformOutput', false);

data = { ...
        rf_05_data, ...
        rf_10_data, ...
        rf_20_data, ...
        rf_40_data, ...
        cmaes_data, ...
        dtscmaes_data
        };
      
datanames = { ...
    'RF 0.05 DTS', ...
    'RF 0.1 DTS', ...
    'RF 0.2 DTS', ...
    'RF 0.4 DTS', ...
    'CMA-ES', ...
    'DTS-CMA-ES' ...
    % 'adaptive DTS-CMA-ES' ...
    };
  
colors = [axisCol; gaussCol; hcCol; pairCol; cmaesCol; dtsCol]/255; 
markers = {rf05Mark, rf10Mark, rf20Mark, rf40Mark, cmaesMark, dtsMark};

plotFuns = [1:24];
plotDims = [2, 3, 5];

clear pdfNames
pdfNames = {};
for d = plotDims
  pdfNames{end+1} = fullfile(plotResultsFolder, sprintf('restrParam_%dD', d));
end

close all
han = relativeFValuesPlot(data, ...
                              'DataNames', datanames, 'DataDims', funcSet.dims, ...
                              'DataFuns', funcSet.BBfunc, 'Colors', colors, ...
                              'PlotFuns', plotFuns, 'PlotDims', plotDims, ...
                              'AggregateDims', false, 'OneFigure', false, ...
                              'Statistic', @median, ...
                              'Quantiles', false(1, 6), ...
                              'AggregateFuns', true, ...
                              'LineSpecification', '-', ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'Markers', markers, ...
                              'PlotGrid', [2, 2], ...
                              'MaxInstances', Inf, ...
                              'FunctionNames', true);
                            
print2pdf(han, pdfNames, 1)

%% Aggregated algorithm comparison: DTS-RF-1, DTS-RF-2 (,DTS-RF-3), CMA-ES, DTS-CMA-ES (, lmm-CMA-ES?) 
% Aggregation of function values in dimensions 2, 3, 5, 10.

data = { ...
        rf_40_gauss_data, ...
        rf_40_res_data, ...
        cmaes_data, ...
        dtscmaes_data, ...
        lmmcmaes_data
        };
      
datanames = { ...
    'SECRET 0.4 DTS', ...
    'SUPPORT 0.4 DTS', ...
    'CMA-ES', ...
    'DTS-CMA-ES', ...
    'lmm-CMA-ES'
    };
  
colors = [rf_40_gaussCol; rf_40_resCol; cmaesCol; dtsCol; lmmCol]/255; 
markers = {gaussMark, resMark, cmaesMark, dtsMark, lmmMark};

plotFuns = [1:24];

% 2D
plotDims = [2];
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
                              'Statistic', 'quantile', ... % @median, ...
                              'Quantiles', true(1, 5), ...
                              'AggregateFuns', false, ...
                              'LineSpecification', '-', ...
                              ... % 'LineWidth', [ 2*ones(1,5), ones(1,7) ], ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'Markers', markers, ...
                              'PlotGrid', [8, 3], ...
                              'ScaleY08', true, ...
                              'MaxInstances', Inf, ...
                              'FunctionNames', true);
                            
print2pdf(han, pdfNames, 1)

% 5D
plotDims = [5];
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
                              'Statistic', 'quantile', ... % @median, ...
                              'Quantiles', true(1, 5), ...
                              'AggregateFuns', false, ...
                              'LineSpecification', '-', ...
                              ... % 'LineWidth', [ 2*ones(1,5), ones(1,7) ], ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'Markers', markers, ...
                              'PlotGrid', [8, 3], ...
                              'ScaleY08', true, ...
                              'MaxInstances', Inf, ...
                              'FunctionNames', true);
                            
print2pdf(han, pdfNames, 1)


%% Algorithm comparison: DTS-RF-1, DTS-RF-2 (,DTS-RF-3), CMA-ES, DTS-CMA-ES (, lmm-CMA-ES?) 
% Aggregation of function values accross dimensions 2, 3, 5, 10.

data = { ...
        rf_40_gauss_data, ...
        rf_40_res_data, ...
        cmaes_data, ...
        dtscmaes_data, ...
        lmmcmaes_data
        };
      
datanames = { ...
    'SECRET 0.4 DTS', ...
    'SUPPORT 0.4 DTS', ...
    'CMA-ES', ...
    'DTS-CMA-ES', ...
    'lmm-CMA-ES'
    };
  
colors = [rf_40_gaussCol; rf_40_resCol; cmaesCol; dtsCol; lmmCol]/255; 
markers = {gaussMark, resMark, cmaesMark, dtsMark, lmmMark};

plotFuns = [1:24];
plotDims = [2, 5];

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
                              'Statistic', @median, ...
                              'Quantiles', false(1, 5), ...
                              'AggregateFuns', true, ...
                              'LineSpecification', '-', ...
                              ... % 'LineWidth', [ 2*ones(1,5), ones(1,7) ], ...
                              'LegendOption', 'split', 'MaxEval', maxEvals, ...
                              'Markers', markers, ...
                              'PlotGrid', [2, 2], ...
                              'ScaleY08', true, ...
                              'MaxInstances', Inf, ...
                              'FunctionNames', true);
                            
print2pdf(han, pdfNames, 1)

%% final clearing
close all
