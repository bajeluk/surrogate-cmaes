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
        [pv, meanRanks] = postHocTest(fValData, 'friedman');
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

      FID = fopen(resultFile, 'w');
      printTableTex(FID, dTable, dims, evaluations, ...
        datanames, pValData, meanRanksData, alpha);
      fclose(FID);
 
      fprintf('Table written to %s\n', resultFile);

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

function printTableTex(FID, table, dims, evaluations, datanames, pVals, ...
  meanRanks, alpha)
% Prints table to file FID

  numOfData = length(datanames);

  nEvals = length(evaluations);
  nDims = length(dims);
  
  % symbol for number of evaluations reaching the best target
  bestSymbol = '\bestFED';
  maxFunEvalsSymbol = '\maxFED';
  maxFunEvalsString = '$250\dm$';
  ftargetString = '10^{-8}';
  NASymbol = '---';

  % representation of evaluation counts as a fraction
  evaluationsString = cell(1, nEvals);
  for e = 1:nEvals
    s = strsplit(strtrim(rats(evaluations(e))), '/');
    if length(s) == 1
      evaluationsString{e} = sprintf('%s', s{1});
    else
      evaluationsString{e} = sprintf('{\\\\LARGE\\\\sfrac{%s}{%s}}', s{1}, s{2});
    end
  end
  
  fprintf(FID, '%%%% for vertically centered cells\n');
  fprintf(FID, '%%\\renewcommand{\\tabularxcolumn}[1]{m{#1}}\n');

  fprintf(FID, '\\begin{table*}[t]\n');
  fprintf(FID, '\\newlength{\\dueltabcolw}\n');
  fprintf(FID, '\\newlength{\\savetabcolsep}\n');
  fprintf(FID, '\\newlength{\\savecmidrulekern}\n');
  fprintf(FID, '\\setlength{\\savetabcolsep}{\\tabcolsep}\n');
  fprintf(FID, '\\setlength{\\savecmidrulekern}{\\cmidrulekern}\n');
  fprintf(FID, '\n');
  fprintf(FID, '\\setlength{\\tabcolsep}{3pt}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{3pt}\n');
  fprintf(FID, '\\setlength{\\dueltabcolw}{\\textwidth-1.85cm-%d\\tabcolsep}\n', 2*(2*numOfData+1));
  fprintf(FID, '\\setlength{\\dueltabcolw}{\\dueltabcolw/%d}\n', 2*numOfData);
  fprintf(FID, '\\centering\n');
  fprintf(FID, '%%\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}X}\n');
  fprintf(FID, '\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}m{\\dueltabcolw}}\n');
  fprintf(FID, '%%\\begin{tabularx}{\\textwidth}{ m{2cm}%s }\n', [repmat('R', 1, 2*numOfData), '']); % 'rr'
  fprintf(FID, '\\begin{tabular}{ m{1.85cm}%s }\n', [repmat('R', 1, 2*numOfData), '']); % 'rr'
  
  for dim = 1:nDims
    if dim == 1
      fprintf(FID, '\\toprule\n');
    else
      fprintf(FID, '\\midrule[\\heavyrulewidth]\n');
    end

    fprintf(FID, '$\\bm{%d\\dm}$ &', dims(dim));

    % header with algorithm names
    % fprintf(FID, strjoin(formatCell('\\\\multicolumn{2}{c}{%s}', datanames), ' & '));
    fprintf(FID, strjoin(formatCell('\\\\multicolumn{2}{c}{\\\\parbox{2\\\\dueltabcolw}{\\\\centering\\\\strut %s\\\\strut}}', datanames), ' & '));
    %fprintf(FID, ' & \\multicolumn{2}{c}{Mean Rank}');
    fprintf(FID, '\\\\\n');

    fprintf(FID, '\\cmidrule(lr){1-1}\n');
    for i = 1:numOfData
      fprintf(FID, '\\cmidrule(lr){%d-%d}\n', 2*i, 2*i+1);
    end

    % header with evaluation numbers
    fprintf(FID, '{\\LARGE\\sfrac{\\nbFEs}{%s}} & ', bestSymbol);
    fprintf(FID, strjoin(repmat(evaluationsString, 1, numOfData), ' & ')); % numOfData + 1
    fprintf(FID, '\\\\\n');
    fprintf(FID, '\\midrule\n');

    for i = 1:numOfData
      % column with algorithm's name
      fprintf(FID, sprintf('%s', datanames{i}));

      for j = 1:numOfData
        for k = 1:nEvals
          dt = table{dim, k};
          pv = pVals{dim, k};

          if i ~= j
            sigdiff = dt(i, j) > dt(j, i) && pv(i, j) < alpha;
            fprintf(FID, ' & ');
            %           if sigdiff
            %             fprintf(FID, '\\hfil');
            %           end
            fprintf(FID, ' %d', dt(i, j));

            if sigdiff
              %             fprintf(FID, '\\newline');
              %             fprintf(FID, ' \\raisebox{2pt}{\\footnotesize\\num{%.2e}}', pv(i, j));
              fprintf(FID, '\\makebox[0pt][l]{$^{\\ast}$}');
            end
          else
            % symbols at diagonal
            fprintf(FID, [' & ', NASymbol]);
          end
        end
      end

      %     % mean rank column
      %     for k = 1:nEvals
      %       mr = meanRanks{k};
      %       fprintf(FID, '& \\multirow{2}{*}{%.2f}', mr(i));
      %     end
      fprintf(FID, '\\\\\n');

    end

    if dim == nDims
      fprintf(FID, '\\bottomrule\n');
    end
  end
  fprintf(FID, '%%\\end{tabularx}\n');
  fprintf(FID, '\\end{tabular}\n');

  % caption printing
  dimString = sprintf('$%d\\\\dm$', dims(1));
  for dim = 2:nDims-1
    dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dims(dim))}, ', ');
  end

  if nDims > 1
    dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dims(nDims))}, ' and ');
  end

  fprintf(FID, ['\\caption{A pairwise comparison of the algorithms in ', dimString, ...
                ' over the BBOB for different evaluation budgets.\n', ... 
                'The number of wins of ith algorithm against jth algorithm ', ...
                'over all benchmark functions is given in ith row and jth column.\n', ...
                'The asterisk marks the better out of each significantly different pair ', ...
                'of algorithms according to the Friedman post-hoc test with the Bergmann-Hommel ', ...
                'correction at significance level $\\alpha=%.2f$.\n', ...
                '%%%s\\ denotes the smallest \\nbFEs\\ at which any of the tested algorithms ', ...
                'reached the target $\\ftarget = %s$.\n', ...
                '%%%s\\ is equal to %s, the overall budget for the experiments.\n', ...
                '}\n'], ...
          alpha, bestSymbol, ftargetString, maxFunEvalsSymbol, maxFunEvalsString ...
  );

  fprintf(FID, '\n');
  fprintf(FID, '\\setlength{\\tabcolsep}{\\savetabcolsep}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{\\savecmidrulekern}\n');

  fprintf(FID, '\\label{tab:duel}\n');
  fprintf(FID, '\\end{table*}\n');

  fprintf(FID, '\\renewcommand{\\tabularxcolumn}[1]{p{#1}}\n');
end

function cellOfStr = formatCell(fmt, cellOfStr)
% map sprintf to each string in a cell array
  cellOfStr = cellfun(@(x) sprintf(fmt, x), cellOfStr, 'UniformOutput', false);
end
