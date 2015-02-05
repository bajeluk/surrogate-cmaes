
exp_id = 'exp_geneEC_06_gplite';

%% Loading filenames

path = ['exp/experiments' filesep exp_id];
mask = [exp_id '_results*.mat'];

algFiles = dir(fullfile(path, mask));
fileIdxs = 1:length(algFiles);

%% Initialization

functions  = [1 2 3 5 6 8 10 11 12 13 14 20 21];
% for 'exp_geneEC_05':
% functions  = [1 2 5 6 8 10 20 21];
dimensions = [2 5 10];

sgParamDef(1).name   = 'evoControlSampleRange';
sgParamDef(1).values = {1, 2};
sgParamDef(2).name   = 'evoControlModelGenerations';
sgParamDef(2).values = {1, 2, 3};

thresholdBounds = { ...
  [0 -8], [0 -8], [0 -8]; ... % f1
  [3 -4], [4  0], [6  1]; ... % f2
  [1  0], [2  0], [2  1]; ... % f3
  [0 -8], [0 -8], [1 -8]; ... % f5
  [1 -6], [2 -4], [3 -3]; ... % f6
  [2 -6], [3 -1], [3  1]; ... % f8
  [3 -3], [4  1], [5  2]; ... % f10
  [3 -3], [3  0], [3  1]; ... % f11
  [3 -2], [4  0], [6  0]; ... % f12
  [1 -3], [2 -1], [2  0]; ... % f13
 [-1 -6], [0 -4], [0 -4]; ... % f14
  [2 -1], [2  0], [2  0]; ... % f20
  [0  0], [1  0], [1  0]  ... % f21
  };

nFiles  = length(fileIdxs);
nParams = length(sgParamDef);
nValues = cell2mat(structMap(sgParamDef, @(x) length(x.values)));
nCombinations = prod(nValues);

evaluations = cell(length(functions), length(dimensions));
wilcoxon    = cell(length(functions), length(dimensions));
cmaesEvaluations = cell(length(functions), length(dimensions));
thresholds  = cell(length(functions), length(dimensions));
speedups    = cell(length(functions), length(dimensions));
filenames   = cell(length(functions), length(dimensions));
params      = cell(length(functions), length(dimensions));

%% Files processing

for fi = fileIdxs
  fname = algFiles(fi).name;
  load(fullfile(path, fname), 'y_evals', 'exp_settings', 'exp_results', 'surrogateParams');

  func   = exp_settings.bbob_function;
  funcId = find(functions == func, 1, 'first');
  dim    = exp_settings.dim;
  dimId  = find(dimensions == dim, 1, 'first');

  % get 10^E thresholds for current func/dim combination
  thisThresholds = thresholdBounds{funcId,dimId}(1):-1:thresholdBounds{funcId,dimId}(2);

  % get the f-evaluations to get to the thresholds
  thresholdsEvals = evaluateYThresholds(y_evals, thisThresholds);

  % calculate (median of ~NaNs) * 1/(% of NaNs)
  expectedEvals = expectedMedian(thresholdsEvals, 2);
  stdEvals      = mapColumns(thresholdsEvals, @std, 2);

  evaluations{funcId, dimId} = [evaluations{funcId, dimId} expectedEvals];
  thresholds{funcId, dimId}  = thisThresholds';

  % save the filenames
  if (isempty(filenames{funcId, dimId}))
    filenames{funcId, dimId} = {fname};
  else
    filenames{funcId, dimId}(end+1) = {fname};
  end
  % save the parameters
  paramVector = zeros(nParams, 1);
  for i = 1:nParams
    paramVector(i) = find(cell2mat(sgParamDef(i).values) == surrogateParams.(sgParamDef(i).name));
  end
  params{funcId, dimId}    = [params{funcId, dimId} paramVector];
  
  % CMA-ES RESULTS
  cmaesFileMask = fullfile(path, 'cmaes_results', [exp_id '_purecmaes_' num2str(func) '_' num2str(dim) 'D_*.mat']);
  cmaesFileEntry = dir(cmaesFileMask);
  load(fullfile(path, 'cmaes_results', cmaesFileEntry.name));
  cmaesEvals = evaluateYThresholds(y_evals, thisThresholds);
  cmaesExpectedEvals = expectedMedian(cmaesEvals, 2);
  cmaesEvaluations{funcId, dimId} = cmaesExpectedEvals;

  % === the final SPEEDUP results ===
  speedups{funcId, dimId} = [speedups{funcId, dimId} cmaesExpectedEvals ./ expectedEvals];

  % === the final Wilcoxon tests ===
  thisWilcoxon = nan(length(thisThresholds),1);
  for th = 1:length(thisThresholds)
    if (sum(~isnan(thresholdsEvals(th,:))) > 2 ...
        && sum(~isnan(cmaesEvals(th,:))) > 2)
      thisWilcoxon(th) = ranksum(thresholdsEvals(th,:), cmaesEvals(th,:), 'tail', 'left');
    end
  end
  wilcoxon{funcId, dimId} = [wilcoxon{funcId, dimId} thisWilcoxon];
end

%% Print the results in LaTeX-friendly mode

% print to stdout
fid = 1;

for fun = 1:length(functions)
  % get the max. number of threshodls in this dimension
  maxThresholds = 0;
  for dim = 1:length(dimensions)
    maxThresholds = max([maxThresholds, length(thresholds{fun, dim})]);
  end
  % lines of text to the table, +1 for the average
  lines = cell(maxThresholds+1, 1);

  fprintf(fid, '% ======  f%d  ======\n', functions(fun));
  % first column spans multiple rows for all the threshodls
  fprintf(fid, '\\hline \\hline\n\\multirow{%d}{*}{\\textbf{f%d}}\n', maxThresholds+1, functions(fun));
  
  for dim = 1:length(dimensions)
    % get the average speedups from each settings
    avgSpeedup = mapColumns(speedups{fun,dim}, @mean);
    % sort the settings according to this avg. speedups
    [~, settingsOrder] = sort(avgSpeedup, 'descend');
    % identify the best and median settings
    bestSettingsIdx = settingsOrder(1);
    mediSettingsIdx  = settingsOrder(floor(length(settingsOrder)/2));

    % add entries to the table
    for th = 1:length(thresholds{fun,dim})
      if (wilcoxon{fun,dim}(th,bestSettingsIdx) <= 0.05)
        speedupFormat = '\\textbf{%.2f}';
      else
        speedupFormat = '%.2f';
      end
      lines{th} = [lines{th} sprintf([' & \\textbf{1e%d} & ' speedupFormat ' & (%.2f) & %.2f   '], thresholds{fun,dim}(th), speedups{fun,dim}(th,bestSettingsIdx), wilcoxon{fun,dim}(th,bestSettingsIdx), speedups{fun,dim}(th,mediSettingsIdx))];
    end
    % generate empty columns for missing thresholds
    for th = (length(thresholds{fun,dim})+1):maxThresholds
      lines{th} = [lines{th} ' & & & &  '];
    end
    lines{end} = [lines{end} sprintf(' & \\textbf{avg:} & %.2f & & %.2f   ', avgSpeedup(bestSettingsIdx), avgSpeedup(mediSettingsIdx))];
  end

  % print the lines
  for l = 1:maxThresholds
    fprintf(fid, '%s \\\\\n', lines{l});
  end
  fprintf(fid, '\\cline{2-13}\n');
  fprintf(fid, '%s \\\\\n', lines{end});
end

