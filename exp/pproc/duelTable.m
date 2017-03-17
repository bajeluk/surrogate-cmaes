function [dTable, ranks] = duelTable(data, varargin)
% [rankTable, ranks] = duelTable(data, settings)
% Creates and prints table containing rankings for different evaluations.
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
%     'Ranking'     - type of ranking (see help createRankingTable)
%                       'tolerant' - equal rank independence
%                       'precise'  - equal ranks shift following ranks
%                       'median'   - equal ranks replaced by medians of
%                                    shifted ranks (from 'precise')
%     'ResultFile'  - file containing resulting table
%     'Statistic'   - statistic of data | string or handle (@mean, @median)
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%     'Alpha'       - significance level for hypothesis testing
%
% Output:
%   rankTable - table of rankings
%   ranks     - rankings for each function and dimension
%
% See Also:
%   createRankingTable, speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  dTable = [];
  if nargin < 1 || isempty(data)
    help duelTable
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
  alpha = defopts(settings, 'alpha', 0.05);
  ftarget = defopts(settings, 'ftarget', 1e-8);
  
  % create ranking table
  extraFields = {'DataNames', 'ResultFile'};
  fieldID = isfield(settings, extraFields);
  createSettings = rmfield(settings, extraFields(fieldID));
  createSettings.Mode = 'target';
  [~, ranks, values] = createRankingTable(data, createSettings);

  nDim = length(dims);
  nEvals = length(evaluations);

  % if there is R-package for computation of p-values
  countPVal = exist('postHocTest', 'file');
  pValData = cell(nDim, nEvals);
  meanRanksData = cell(nDim, nEvals);
  dTable = cell(nDim, nEvals);

  for d = 1:nDim
    for e = 1:nEvals
      if countPVal
        fValData = cell2mat(arrayfun(@(x) values{x, d}(e, :), BBfunc, 'UniformOutput', false)');
        [pv, meanRanks] = postHocTest(fValData, 'friedman', sprintf('/tmp/friedman_%d_%d.out.csv', d, e));
        pValData{d, e} = pv;
        meanRanksData{d, e} = meanRanks;
      end
      rankData = cell2mat(arrayfun(@(x) ranks{x, d}(e, :), BBfunc, 'UniformOutput', false)');
      dTable{d, e} = createDuelTable(rankData);
    end
  end
  
  % print table
  switch tableFormat
    % prints table to latex file
    case {'tex', 'latex'}
      if ~exist(resultFolder, 'dir')
        mkdir(resultFolder)
      end

      if nDim > 1
        resultFile = arrayfun(@(x) [resultFile(1:end-4), '_', num2str(x), ...
          'D', resultFile(end-3:end)], dims, 'UniformOutput', false);
      else
        resultFile{1} = resultFile;
      end

      for d = 1:nDim
        FID = fopen(resultFile{d}, 'w');
        printTableTex(FID, dTable(d, :), dims(d), evaluations, datanames, pValData(d, :), ...
          meanRanksData(d, :), alpha, ftarget);
        fclose(FID);
        
        fprintf('Table written to %s\n', resultFile{d});
      end

    otherwise
      error('Format ''%s'' is not supported.', tableFormat)
  end

end

function dt = createDuelTable(ranks)
% create duel table
% rank of data in row is lower than rank of data in column
  [nFun, nData] = size(ranks);

  dt = zeros(nData);
  for f = 1:nFun
    for dat = 1:nData
      id = ranks(f, dat) < ranks(f, :);
      dt(dat, :) = dt(dat, :) + id;
    end
  end
end

function printTableTex(FID, table, dim, evaluations, datanames, pVals, meanRanks, alpha, ftarget)
% Prints table to file FID
  FID =1;

  numOfData = length(datanames);

  nEvals = length(evaluations);
  
  % symbol for number of evaluations reaching the best target
  bestSymbol = '\\bestFED';
  maxFunEvalsString = '\\maxfunevals';

  fprintf(FID, '\\begin{table}\n');
  fprintf(FID, '\\centering\n');
  fprintf(FID, '\\begin{tabular}{ l%s }\n', [repmat('c', 1, 2*numOfData), 'rr']);
  fprintf(FID, '\\toprule');
  fprintf(FID, '\\textbf{%dD} &', dim);
  %fprintf(FID, '\\multicolumn{%d}{c}{Score} & Mean ranks\\\n', numOfData);
  %fprintf(FID, '\\midrule\n');

  % make datanames equally long
  datanames = sameLength(datanames);

  % header with algorithm names
  fprintf(FID, strjoin(formatCell('\\\\multicolumn{2}{c}{%s}', datanames), ' & '));
  fprintf(FID, ' Mean Rank &\\\\\n');
  fprintf(FID, '\\midrule\n');

  % header with evaluation numbers
  fprintf(FID, '# FE');
  fprintf(FID, ' & %d', evaluations);
  fprintf(FID, '\\\\\n');
  fprintf(FID, '\\midrule\n');

  for i = 1:numOfData
    % column with algorithm's name
    fprintf(FID, '\\multirow{2}{*}{%s}', datanames{i});
    
    % a row with mutual score
    for j = 1:numOfData
      for k = 1:nEvals
        dt = table{k};
        fprintf(FID, ' & %d', dt(i, j));
      end
    end
    fprintf(FID, '\\\\\n');

    
    % mean rank column
    for k = 1:nEvals
      mr = meanRanks{k};
      fprintf(FID, '& \\multirow{2}{*}{%.2f}', mr(i));
    end
    fprintf(FID, '\\\\\n');

    % a subrow with pvalues
    fprintf(FID, ' ');
    for j = 1:numOfData
      for k = 1:nEvals
        dt = table{k};
        pv = pVals{k};
        if dt(i, j) > dt(j, i) && pv(i, j) < alpha
          fprintf(FID, ' & {\\footnotesize %.2e}', pv(i, j));
        else
          fprintf(FID, ' & {}');
        end
      end
    end
  end
  fprintf(FID, '\\bottomrule\n');
  fprintf(FID, '\\end{tabular}\n');

  % caption printing
  % fprintf(FID, '\\vspace{1mm}\n');
  fprintf(FID, ['\\caption{Multi-comparison of algorithms in %dD. \n', ...
                '%s denotes the smallest FE/D at which any ', ...
                'of the tested algorithms reached the target \\Delta_f^\\text{med} = %.0e ', ...
                'or %s if the target ', ...
                'was not reached by any algorithms. \n', ...
                'The number of wins of ith algorithm over jth algorithm ', ...
                'for different FE/D ', ...
                'over all benchmark functions are given in ith row and jth column. \n', ...
                'P-value of Friedman post-hoc test with significance level \\alpha=%.2f ', ...
                'is given for the winning algorithm in each significantly different pair.\n'], ...
                dim, bestSymbol, ftarget, maxFunEvalsString, alpha);

  fprintf(FID, '\\label{tab:duel%d}\n', dim);
  fprintf(FID, '\\end{table}\n');
  
end

function cellOfStr = sameLength(cellOfStr)
% returns cell array of strings with added empty space to the same length
  maxLength = max(cellfun(@length, cellOfStr));
  cellOfStr = cellfun(@(x) [x, repmat(' ', 1, maxLength - length(x))], cellOfStr, 'UniformOutput', false);
end

function cellOfStr = formatCell(fmt, cellOfStr)
% map sprintf to each string in a cell array
  cellOfStr = cellfun(@(x) sprintf(fmt, x), cellOfStr, 'UniformOutput', false);
end
