function [table, ranks] = rankingTable(data, varargin)
% [table, ranks] = rankingTable(data, settings)
% Creates table containing rankings for different evaluations.
%
% Input:
%   data      - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataNames'   - cell array of data names (e.g. names of algorithms)
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%     'Evaluations' - evaluations chosen to count
%     'Statistic'   - statistic of data | string or handle (@mean, 
%                       @median)
%     'ResultFile'  - file containing resulting table
%
% Output:
%   table - table of rankings
%   ranks - rankings for each function and dimension
%
% See Also:
%   speedUpPlot, speedUpPlotCompare, dataReady

  % initialization
  table = [];
  if nargin < 1 || isempty(data)
    help relativeFValuesPlot
    return
  end
  if isstruct(varargin)
    settings = varargin;
  else
    % keep cells as cells due to struct command
    vararCellId = cellfun(@iscell, varargin);
    varargin(vararCellId) = {varargin(vararCellId)};
    settings = struct(varargin{:});
  end
  numOfData = length(data);
  datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  dims    = defopts(settings, 'TableDims', funcSet.dims);
  BBfunc  = defopts(settings, 'TableFuns', funcSet.BBfunc);
  evaluations = defopts(settings, 'Evaluations', [20 40 80]);
  statistic = defopts(settings, 'Statistic', @median);
  defResultFolder = fullfile('exp', 'pproc', 'tex');
  resultFile = defopts(settings, 'ResultFile', fullfile(defResultFolder, 'rankTable.tex'));
  fileID = strfind(resultFile, filesep);
  resultFolder = resultFile(1 : fileID(end) - 1);
  if ischar(statistic)
    if strcmp(statistic, 'quantile')
      statistic = @(x, dim) quantile(x, [0.25, 0.5, 0.75], dim);
    else
      statistic = str2func(statistic);
    end
  end

  % get function and dimension IDs
  dimInvIds = inverseIndex(funcSet.dims);
  dimIds = dimInvIds(dims);
  funcInvIds = inverseIndex(funcSet.BBfunc);
  funcIds = funcInvIds(BBfunc);

  if ~all(dimIds)
    fprintf('Wrong dimesion request!\n')
  end
  if ~all(funcIds)
    fprintf('Wrong function request!\n')
  end

  % count means
  useMaxInstances = 15;
  data_stats = cellfun(@(D) gainStatistic(D, dimIds, funcIds, ...
                            'MaxInstances', useMaxInstances, ...
                            'AverageDims', false, ...
                            'Statistic', statistic, ...
                            'SuppWarning', true), ...
                       data, 'UniformOutput', false);
                     
  % gain rankings
  nFunc = length(funcIds);
  nDims = length(dimIds);
  nEvals = length(evaluations);
  ranks = cell(nFunc, nDims);
  for f = 1:nFunc
    for d = 1:nDims
      notEmptyData = inverseIndex(arrayfun(@(x) ~isempty(data_stats{x}{f,d}), 1:numOfData));
      for e = 1:nEvals
        thisData = cell2mat(arrayfun(@(x) data_stats{x}{f,d}(evaluations(e)), notEmptyData, 'UniformOutput', false));
        thisData = max(thisData, 1e-8 * ones(size(thisData)));
        [~, ~, I] = unique(thisData);
%         [~, I2] = sort(cell2mat(arrayfun(@(x) data_stats{x}{f,d}(e), notEmptyData, 'UniformOutput', false)));
        ranks{f,d}(e, notEmptyData) = I';
      end
    end
  end

  % aggregate ranks accross functions
  for dat = 1:numOfData
    for d = 1:nDims
      for r = 1:numOfData
        for e = 1:nEvals
          table{dat, d}(r, e) = sum(arrayfun(@(x) ranks{x,d}(e, dat) == r, 1:nFunc));
        end
      end
    end
  end
  
  % print table
  if ~exist(resultFolder, 'dir')
    mkdir(resultFolder)
  end
  FID = fopen(resultFile, 'w');
  fprintf(FID, '\\begin{table}\n');
  fprintf(FID, '\\centering\n');
  fprintf(FID, '\\begin{tabular}[pos]{| l %s |}\n', repmat('|c',1, nEvals*(nDims+1)));
  fprintf(FID, '\\hline\n');
  fprintf(FID, '{} ');
  for d = 1:nDims
    fprintf(FID, '& \\multicolumn{%d}{|c|}{%dD} ', nEvals, dims(d));
  end
  fprintf(FID, '& \\multicolumn{%d}{|c|}{$\\sum$} \\\\\n', nEvals);
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
  % data rows
  for dat = 1:numOfData
    printString = '';
    % dimensions
    for d = 1:nDims
      for e = 1:nEvals
        printString = [printString, ' & ', num2str(table{dat,d}(1,e))];
      end
    end
    % sum
    for e = 1:nEvals
      printString = [printString, ' & ', num2str(sum(arrayfun(@(x) table{dat,x}(1,e), 1:nDims)))];
    end
    fprintf(FID, '%s%s \\\\\n', datanames{dat}, printString);
    fprintf(FID, '\\hline\n');
  end
  fprintf(FID, '\\end{tabular}\n');
  printString = num2str(evaluations(1));
  for e = 2:nEvals
    printString = [printString, ', ', num2str(evaluations(e))];
  end
  fprintf(FID, '\\vspace{1mm}\n');
  fprintf(FID, '\\caption{Counts of the 1st ranks according to the lowest achieved $\\Delta_f$ for different FE/D = \\{%s\\}.}\n', printString);
  fprintf(FID, '\\label{tab:fed}\n');
  fprintf(FID, '\\end{table}\n');

  fclose(FID);
  
  fprintf('Table written to %s\n', resultFile);

end

function cellOfStr = sameLength(cellOfStr)
% returns cell array of strings with added empty space to the same length
  maxLength = max(cellfun(@length, cellOfStr));
  cellOfStr = cellfun(@(x) [x, repmat(' ', 1, maxLength - length(x))], cellOfStr, 'UniformOutput', false);
end