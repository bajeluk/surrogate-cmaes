%% EC 2017 journal CPU Timing table

%% load data

% checkout file containing all loaded data
tmpFName = fullfile('/tmp', 'ec2017timing.mat');
if (exist(tmpFName', 'file'))
  load(tmpFName);
else

% needed function and dimension settings
funcSet.BBfunc = 1:24;
funcSet.dims = [2, 3, 5, 10, 20];
maxEvals = 250;

% folder for results
actualFolder = pwd;
articleFolder = fullfile(actualFolder(1:end - 1 - length('src/surrogate-cmaes')), 'Documents', 'latex_scmaes', 'ec2016paper');
plotResultsFolder = fullfile(articleFolder, 'img');
tableFolder = fullfile(articleFolder, 'data');
[~,~] = mkdir(plotResultsFolder);
[~,~] = mkdir(tableFolder);

% path settings
exppath = fullfile('exp', 'experiments');

cmaes_2pop_path = fullfile(exppath, 'exp_CMAES_2pop');
scmaes_10D_path = fullfile(exppath, 'exp_geneEC_10');
scmaes_20D_path = fullfile(exppath, 'exp_geneEC_10_20D');
% TODO: change to exp_doubleEC_26_1model
dts005_path = fullfile(exppath, 'exp_doubleEC_26');
% TODO: change to exp_doubleEC_26_1model_adapt04
%dts4crit_path = fullfile(exppath, 'exp_doubleEC_26_4crit');
dts_adapt_path = fullfile(exppath, 'exp_doubleEC_26_adapt04');
maes_10D_path = fullfile(exppath, 'exp_maesEC_16_2_10_cmaes_10D');
maes_20D_path = fullfile(exppath, 'exp_maesEC_16_2_10_cmaes_20D');
gpop_10D_path = fullfile(exppath, 'exp_gpop_16_10D');
gpop_20D_path = fullfile(exppath, 'exp_gpop_16_20D');
bobyqa_path = fullfile(exppath, 'exp_BOBYQA_03');
fmincon_pure_path = fullfile(exppath, 'exp_fmincon_03');

% load data
dataFolders = {cmaes_2pop_path; ...
               scmaes_10D_path; ...
               scmaes_20D_path; ...
               dts005_path; ...
               dts_adapt_path; ...
               maes_path; ...
               gpop_10D_path; ...
               gpop_20D_path; ...
               bobyqa_path; ...
               fmincon_pure_path ...
               };

[evals, settings, exp_results] = catEvalSet(dataFolders, funcSet);

% find ids in settings
clear findSet
findSet.evoControl = 'none';
cmaes_2pop_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'generation';
findSet.modelType = 'gp';
findSet.evoControlModelGenerations = 5;
scmaes_gp_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'generation';
findSet.modelType = 'rf';
findSet.evoControlModelGenerations = 5;
scmaes_rf_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'generation';
findSet.modelType = 'gp';
findSet.evoControlModelGenerations = 1;
scmaes_1_gp_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'generation';
findSet.modelType = 'rf';
findSet.evoControlModelGenerations = 1;
scmaes_1_rf_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'doubletrained';
findSet.DTAdaptive_updateRateDown = 'obj.updateRate';
dts_adapt_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'doubletrained';
findSet.DTAdaptive_updateRate = 0.2;
findSet.DTAdaptive_updateRateDown = 'obj.updateRate';
dts_adapt_2_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'doubletrained';
findSet.DTAdaptive_updateRate = 0.4;
findSet.DTAdaptive_updateRateDown = 'obj.updateRate';
dts_adapt_4_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'doubletrained';
findSet.modelOpts.predictionType = 'fvalues';
dts_critfvalues_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'doubletrained';
findSet.modelOpts.predictionType = 'sd2';
dts_critsd2_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'doubletrained';
findSet.modelOpts.predictionType = 'ei';
dts_critei_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'doubletrained';
findSet.modelOpts.predictionType = 'expectedRank';
dts_criterde_Id = getStructIndex(settings, findSet);

clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'fvalues';
maes_Id = getStructIndex(settings, findSet);
clear findSet
findSet.evoControl = 'maes';
findSet.modelOpts.predictionType = 'poi';
maes_poi_Id = getStructIndex(settings, findSet);

clear findSet
findSet.gpop_funHistLen = 5;
gpop_Id = getStructIndex(settings, findSet);
clear findSet
findSet.Algorithm = 'interior-point';
fmincon_pure_Id = getStructIndex(settings, findSet);
clear findSet
findSet.rho_beg = 3;
bobyqa_Id = getStructIndex(settings, findSet);

all_Ids = [cmaes_2pop_Id, dts_critfvalues_Id, dts_critsd2_Id, dts_critei_Id, dts_criterde_Id, ...
  dts_adapt_Id, dts_adapt_2_Id, dts_adapt_4_Id, ...
  scmaes_gp_Id, scmaes_rf_Id, scmaes_1_gp_Id, scmaes_1_rf_Id, ...
  maes_Id, maes_poi_Id, ...
  gpop_Id, fmincon_pure_Id, bobyqa_Id];

dts005_Id = setdiff(1:length(settings), all_Ids);

% extract data
scmaes_gp_data = evals(:, :, scmaes_gp_Id);
dts005_data = evals(:, :, dts005_Id);
maes_data = evals(:, :, maes_Id);
gpop_data = evals(:, :, gpop_Id);
bobyqa_data = evals(:, :, bobyqa_Id);
fmincon_pure_data = evals(:, :, fmincon_pure_Id);

% extract experiment info such as timing
cmaes_2pop_info = exp_results(:, :, cmaes_2pop_Id);
scmaes_gp_info = exp_results(:, :, scmaes_gp_Id);
dts005_info = exp_results(:, :, dts005_Id);
dts_adapt_info = exp_results(:, :, dts_adapt_Id);
maes_info = exp_results(:, :, maes_Id);
gpop_info = exp_results(:, :, gpop_Id);
bobyqa_info = exp_results(:, :, bobyqa_Id);
fmincon_pure_info = exp_results(:, :, fmincon_pure_Id);

if (~exist(tmpFName, 'file'))
  save(tmpFName);
end

end

%% CPU timing table

info_all = {...
  scmaes_gp_info, ...
  dts005_info, ...
  dts_adapt_info, ...
  maes_info, ...
  gpop_info, ...
  cmaes_2pop_info, ...
  bobyqa_info, ...
  fmincon_pure_info ...
  };

datanames = {...
  'S-CMA-ES', ...
  '$0.05/2\mathrm{pop}$ DTS-CMA-ES', ...
  'adaptive DTS-CMA-ES', ...
  'MA-ES', ...
  'GPOP', ...
  'CMA-ES $2\mathrm{pop}$', ...
  'BOBYQA', ...
  'fmincon' ...
  };

timing = cpuTimingTable(info_all, ...
                   'DataNames', datanames, ...
                   'DataFuns', funcSet.BBfunc, 'DataDims', funcSet.dims);
tab_data = cell(length(info_all), length(funcSet.dims)+1);

for alg = 1:length(info_all)
  tab_data{alg, 1} = datanames{alg};
  tab_data(alg, 2:end) = num2cell(mean(timing(:, :, alg)));
end

header = arrayfun(@(x)sprintf('D%d', x), funcSet.dims, 'UniformOutput', false);
header = [{'Algorithm'} header];

timingTable = cell2table(tab_data, 'VariableNames', header);

lt = LatexTable(timingTable);
lt.headerRow = arrayfun(@(x)sprintf('{%d-D}', x), funcSet.dims, 'UniformOutput', false);
lt.headerRow = [{'Algorithm'} lt.headerRow];
lt.opts.tableColumnAlignment = num2cell(['l' ...
  repmat('S[table-number-alignment=center,table-sign-mantissa,table-figures-integer=1,table-figures-decimal=1,table-figures-exponent=1]', ...
  1, length(funcSet.dims))]);
lt.opts.numericFormat = '%2.1e';
lt.opts.booktabs = 1;
lt.opts.latexHeader = 0;

% print the table
fid = fopen(fullfile(tableFolder, 'timing.tex'), 'w');
latexRows = lt.toStringRows(lt.toStringTable());
for i = 1:length(latexRows)
  fprintf(fid, '%s\n', latexRows{i});
end
fclose(fid);