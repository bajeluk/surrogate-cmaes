function [rankTable, ranks] = rankingTable(data, varargin)
% [rankTable, ranks] = rankingTable(data, settings)
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
%     'Target'      - target value to reach
%
% Output:
%   rankTable - table of rankings
%   ranks     - rankings for each function and dimension
%
% See Also:
%   createRankingTable, speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  rankTable = [];
  if nargin < 1 || isempty(data)
    help rankingTable
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
  tableRank = defopts(settings, 'Rank', 1);
  targetValue = defopts(settings, 'Target', 1e-8);
  evaluations = defopts(settings, 'Evaluations', [20 40 80]);
  defResultFolder = fullfile('exp', 'pproc', 'tex');
  resultFile = defopts(settings, 'ResultFile', fullfile(defResultFolder, 'rankTable.tex'));
  fileID = strfind(resultFile, filesep);
  resultFolder = resultFile(1 : fileID(end) - 1);
  
  % create ranking table
  extraFields = {'DataNames', 'ResultFile'};
  fieldID = isfield(settings, extraFields);
  createSettings = rmfield(settings, extraFields(fieldID));
  % createSettings.Mode = 'target';
  [rankTable, ranks] = createRankingTable(data, createSettings);
  
  % print table
  switch tableFormat
    % prints table to figure
    case 'figure'
      nEvals = length(evaluations);
      nDims = length(dims);
      maxLengthData = max(cellfun(@length, datanames));
      
      evalRow = repmat(evaluations, [1, length(dims)+1]);      
      % numbers should have normalized format
      % choose how to display them:
      % a) transform to text
      % table = arrayfun(@(x) sprintf('%g', x), table, 'UniformOutput', false);
      % maxLengthNumber = max(max(cellfun(@length, table)));
      % b) round values
      if tableRank == 1
        rankTable = round(rankTable);
      else
      % c) multiply by ten
        rankTable = 10*(rankTable);
      end
      publicTable = [evalRow; rankTable];
      maxLengthNumber = max(max(arrayfun(@(x) ceil(log10(x)), publicTable)));
      % column width is number-length dependent
      colWidth = 20 + 5*maxLengthNumber;
      tableSize = [11*(2+maxLengthData) + nEvals*(nDims+1)*colWidth, 20*(numOfData+2)];
      
      colBase = arrayfun(@(x) repmat({[num2str(x), 'D']}, [1, nEvals]), dims, 'UniformOutput', false);
      colName = [[colBase{:}], repmat({'SUM'}, [1, nEvals])];
      rowName = [{'FE/D'}, datanames];
      f = figure('Position', [0, 0, tableSize]);
      rankTable = uitable(f, 'Data', publicTable, ...
                 'ColumnName', colName, ...
                 'RowName', rowName, ...
                 'ColumnWidth', {colWidth}, ...
                 'Position', [1, 1, tableSize]);
      
    % prints table to latex file
    case {'tex', 'latex'}
      if ~exist(resultFolder, 'dir')
        mkdir(resultFolder)
      end
      FID = fopen(resultFile, 'w');
      printTableTex(FID, rankTable, dims, evaluations, datanames, length(BBfunc), targetValue)
      fclose(FID);

      fprintf('Table written to %s\n', resultFile);
  end

end

function printTableTex(FID, table, dims, evaluations, datanames, nFunc, targetValue)
% Prints table to file FID

  [numOfData, nColumns] = size(table);
  nDims = length(dims);
  nEvals = length(evaluations);
  
  fprintf(FID, '\\begin{table}\n');
  fprintf(FID, '\\centering\n');
  fprintf(FID, '\\begin{tabular}[pos]{ l %s }\n', repmat([' |', repmat(' c', 1, nEvals)], 1, nDims+1));
  fprintf(FID, '\\hline\n');
  fprintf(FID, '{} ');
  for d = 1:nDims
    fprintf(FID, '& \\multicolumn{%d}{c|}{%dD} ', nEvals, dims(d));
  end
  fprintf(FID, '& \\multicolumn{%d}{c}{$\\sum$} \\\\\n', nEvals);
  printString = '';
  for d = 1:nDims + 1
    for e = 1:nEvals
      printString = [printString, ' & ', num2str(evaluations(e))];
    end
  end
  fprintf(FID, 'FE/D %s \\\\\n', printString);
  fprintf(FID, '\\hline\n');
  % make datanames equally long
  datanames = sameLength(datanames);
  % find max sums of ranks
  maxTableRanks = max(table);
  minTableRanks = min(table);
  % data rows
  for dat = 1:numOfData
    printString = '';
    % columns
    for col = 1:nColumns
      sumRank = table(dat, col);
      if sumRank == minTableRanks(col)
        % print best data in bold
        printString = [printString, ' & ', '\textbf{', num2str(sumRank), '}'];
      else
        printString = [printString, ' & ', num2str(sumRank)];
      end
    end
    fprintf(FID, '%s%s \\\\\n', datanames{dat}, printString);
  end
  fprintf(FID, '\\hline\n');
  fprintf(FID, '\\end{tabular}\n');
  % evaluation numbers
  evalString = arrayString(evaluations, ',');
  % dimension numbers 
  dimString = arrayString(dims, ',');
  % caption printing
  fprintf(FID, '\\vspace{1mm}\n');
  fprintf(FID, ['\\caption{Counts of the 1st ranks from %d benchmark functions according to the lowest achieved ', ...
                '$\\Delta_f^\\text{med}$ for different FE/D = \\{%s\\} and dimensions D = \\{%s\\}. ', ...
                'Ties of the 1st ranks are counted for all respective algorithms. ', ...
                'The ties often occure when $\\Delta f_T = %s$ is reached (mostly on f1 and f5).}\n'], ...
                nFunc, evalString, dimString, num2tex(targetValue));
               
  fprintf(FID, '\\label{tab:fed}\n');
  fprintf(FID, '\\end{table}\n');
  
end

function str = arrayString(vector, delimiter)
% returns string containing 'vector' elements separated by 'delimiter'
  str = num2str(vector(1));
  for e = 2:length(vector);
    str = [str, delimiter, ' ', num2str(vector(e))];
  end
end

function cellOfStr = sameLength(cellOfStr)
% returns cell array of strings with added empty space to the same length
  maxLength = max(cellfun(@length, cellOfStr));
  cellOfStr = cellfun(@(x) [x, repmat(' ', 1, maxLength - length(x))], cellOfStr, 'UniformOutput', false);
end
