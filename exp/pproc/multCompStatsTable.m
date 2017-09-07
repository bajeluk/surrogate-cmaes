function [stats, meanRanks] = multCompStatsTable(data, varargin)
% statsTable = multCompStatsTable
% Creates and prints table with multicomparison statistics, such
%   as F_F values and mean ranks from Friedman or Iman-Davenport test.
%
% Input:
%   data     - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataNames'   - cell array of data names (e.g. names of algorithms)
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'Evaluations' - evaluations chosen to count
%     'Format'      - table format | ('tex', 'figure')
%     'ResultFile'  - file containing resulting table
%     'Statistic'   - statistic of data | string or handle (@mean, @median)
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%
% Output:
%   stats     - Friedman statistics for each dimension and each
%               evaluations number
%   meanRanks - mean ranks over all functions for each dimension and each
%               evaluations number

  % initialization
  stats = [];
  if nargin < 1 || isempty(data)
    help multCompStatsTable
    return
  end
  settings = settings2struct(varargin);

  numOfData = length(data);
  datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  tableFormat = defopts(settings, 'Format', 'tex');
  dims    = defopts(settings, 'TableDims', funcSet.dims);
  BBfunc  = defopts(settings, 'TableFuns', funcSet.BBfunc);
  evaluations = defopts(settings, 'Evaluations', [1/3, 1]);
  defResultFolder = fullfile('exp', 'pproc', 'tex');
  resultFile = defopts(settings, 'ResultFile', fullfile(defResultFolder, 'duelTable.tex'));
  fileID = strfind(resultFile, filesep);
  resultFolder = resultFile(1 : fileID(end) - 1);
  printHeader = defopts(settings, 'PrintHeader', true);
  
  nDim = length(dims);
  nEvals = length(evaluations);

  meanRanks = cell(nDim, nEvals);
  stats = cell(nDim, nEvals);
  pvals = cell(nDim, nEvals);

  % create tables with values at different budgets of evaluations
  extraFields = {'DataNames', 'ResultFile'};
  fieldID = isfield(settings, extraFields);
  createSettings = rmfield(settings, extraFields(fieldID));
  createSettings.Mode = 'target';
  [~, ~, values] = createRankingTable(data, createSettings);

  for d = 1:nDim
    for e = 1:nEvals
      fValData = cell2mat(arrayfun(@(x) values{x, d}(e, :), BBfunc, 'UniformOutput', false)');
      [~, mr] = postHocTest(fValData, 'friedman', 'shaffer');
      [p, stat] = multipleComparisonTest(fValData, 'iman');
      meanRanks{d, e} = mr;
      stats{d, e} = stat;
      pvals{d, e} = p;
    end
  end

  % print table
  switch tableFormat
    % prints table to latex file
    case {'tex', 'latex'}
      if ~exist(resultFolder, 'dir')
        mkdir(resultFolder);
      end
      FID = fopen(resultFile, 'w');
      printTableTex(FID, stats, pvals, '$F_F$', meanRanks, dims, evaluations, ...
        datanames, printHeader);
      fclose(FID);

      fprintf('Table written to %s\n', resultFile);
    otherwise
      error('Format ''%s'' is not supported.', tableFormat);
  end
end

function printTableTex(FID, stats, pvals, statsSymbol, meanRanks, dims, evaluations, ...
  datanames, printHeader)

  numOfData = length(datanames);
  nDims = length(dims);
  nEvals = length(evaluations);

  % symbol for number of evaluations reaching the best target
  bestSymbol = '\bestFED';
  maxFunEvalsSymbol = '\maxFED';
  
  % representation of evaluation counts as a fraction
  evaluationsString = cell(1, nEvals);
  for e = 1:nEvals
    s = strsplit(strtrim(rats(evaluations(e))), '/');
    if length(s) == 1
      evaluationsString{e} = s{1};
    else
      evaluationsString{e} = sprintf('{\\\\large\\\\sfrac{%s}{%s}}', s{1}, s{2});
    end
  end

  if printHeader
    fprintf(FID, '\\begin{table}[t]\n');
    fprintf(FID, '\\centering\n');
  end

  fprintf(FID, '\\begin{tabular}{ l%s }\n', repmat('r', 1, nDims * nEvals));
  fprintf(FID, '\\toprule\n');

  % dimensionality header
  fprintf(FID, 'Dim');
  fprintf(FID, sprintf('& \\\\multicolumn{2}{c}{$%d\\\\dm$}', dims));
  fprintf(FID, '\\\\\n');

  fprintf(FID, '\\cmidrule(lr){1-1}\n');
  for i = 1:nDims
    fprintf(FID, '\\cmidrule(lr){%d-%d}\n', 2*i, 2*i+1);
  end

  % evaluations header
  fprintf(FID, '{\\LARGE\\sfrac{\\nbFEs}{%s}} & ', bestSymbol);
  fprintf(FID, strjoin(repmat(evaluationsString, 1, nDims), ' & '));
  fprintf(FID, '\\\\\n');
  
  fprintf(FID, '\\midrule\n');

  mins = cell(nDims, nEvals);
  mins(1:nDims, 1:nEvals) = {Inf};

  for d = 1:nDims
    for e = 1:nEvals
      r = meanRanks{d, e};
      [~, ind] = min(r);
      mins{d, e} = ind;
    end
  end
  
  % rows
  for i = 1:numOfData
    % the column with algorithm's name
    fprintf(FID, sprintf('%s', datanames{i}));
    
    for d = 1:nDims
      for e = 1:nEvals
        mr = meanRanks{d, e};
        if i == mins{d, e}
          fprintf(FID, ' & $\\bm{%.2f}$', mr(i));
        else
          fprintf(FID, ' & $%.2f$', mr(i));
        end
      end
    end
    fprintf(FID, '\\\\\n');
  end
  
  fprintf(FID, '\\cmidrule(lr){1-1}\n');
  for i = 1:nDims
    fprintf(FID, '\\cmidrule(lr){%d-%d}\n', 2*i, 2*i+1);
  end

  % additional row with the statistics
  fprintf(FID, statsSymbol);
  for d = 1:nDims
    for e = 1:nEvals
      fprintf(FID, ' & $%.2f$', stats{d, e});
    end
  end

  % and with the p-values
  fprintf(FID, '\\\\\np-val');
  for d = 1:nDims
    for e = 1:nEvals
      fprintf(FID, ' & $%.2e$', pvals{d, e});
    end
  end

  fprintf(FID, '\\\\\n');
  fprintf(FID, '\\bottomrule\n');
  fprintf(FID, '\\end{tabular}\n');
  
  if printHeader
    % caption
    fprintf(FID, ['\\caption{Mean rank of each algorithm over the BBOB ', ...
      'and the Iman-Davenport statistic ', ...
      'for the %d considered combinations of ', ...
      'dimensionalities and evaluation budgets.}\n'], nEvals * nDims);
    fprintf(FID, '\\label{tab:stat}\n');
    fprintf(FID, '\\end{table}\n');
  end
end

