function prtSignifTable(data, varargin)
% significanceTable(data, varargin)

  if nargin < 1 || isempty(data)
    help prtSignifTable
    return
  end
  settings = settings2struct(varargin);

  par.rowNames = defopts(settings, 'RowNames', {});
  par.rowValName = defopts(settings, 'RowValName', '');
  par.rowGroups = defopts(settings, 'RowGroups', {});
  par.rowGroupNum = defopts(settings, 'RowGroupNum', []);
  par.colNames = defopts(settings, 'ColNames', {});
  par.colValName = defopts(settings, 'ColValName', '');
  par.colGroups = defopts(settings, 'ColGroups', {});
  par.colGroupNum = defopts(settings, 'ColGroupNum', []);
  defResultFolder = fullfile('exp', 'pproc', 'tex');
  resultFile = defopts(settings, 'ResultFile', ...
                       fullfile(defResultFolder, 'signifTable.tex'));
  [resultFolder, par.resultLabel] = fileparts(resultFile); 
  par.alpha = defopts(settings, 'Alpha', 0.05);
  par.printHeader = defopts(settings, 'PrintHeader', true);
  par.tableWidth = defopts(settings, 'TableWidth', '\textwidth');
  par.oneColumn = defopts(settings, 'OneColumn', false);
  par.caption = defopts(settings, 'Caption', 'Caption text.');
  
  % check values
  assert(numel(par.rowNames) == size(data, 1) || numel(par.rowNames) == 0, ...
    'Number of row names and number of rows are not identical')
  assert(numel(par.colNames) == size(data, 2) || numel(par.colNames) == 0, ...
    'Number of column names and number of columns are not identical')
  assert(numel(par.rowGroups) == numel(par.rowGroupNum), ...
    'Number of row groups and number of group rows are not identical')
  assert(numel(par.colGroups) == numel(par.colGroupNum), ...
    'Number of column groups and number of group columns are not identical')
  
  % print table
  [~, ~] = mkdir(resultFolder);
  FID = fopen(resultFile, 'w');
  printTableTex(FID, data, par);
  fclose(FID);
 
  fprintf('Table written to %s\n', resultFile);
                     
end

function printTableTex(FID, data, par)
% Prints table to file FID

  [nRows, nCols] = size(data);
  nRowGroups = numel(par.rowGroups);
  nColGroups = numel(par.colGroups);
  
  % constants
  rowShrink = 2;
  
  % header printing
  if par.printHeader
    if par.oneColumn
      fprintf(FID, '\\begin{table}[t]\n');
    else
      fprintf(FID, '\\begin{table*}[t]\n');
    end
    
    % caption printing
    fprintf(FID, '\\caption{');
    fprintf(FID, '%s', par.caption);
    fprintf(FID, '}\n');
  end

  % define lengths
  fprintNewLength(FID, 'tabcolw');
  fprintNewLength(FID, 'savetabcolsep')
  fprintNewLength(FID, 'savecmidrulekern');
  fprintNewLength(FID, 'headcolw');
  fprintNewLength(FID, 'groupheadcolw');
  fprintf(FID, '\\setlength{\\savetabcolsep}{\\tabcolsep}\n');
  fprintf(FID, '\\setlength{\\savecmidrulekern}{\\cmidrulekern}\n');
  fprintf(FID, '\n');
  fprintf(FID, '\\setlength{\\headcolw}{1.5cm}\n');
  fprintf(FID, '\\setlength{\\groupheadcolw}{0.8cm}\n');
  fprintf(FID, '\\setlength{\\tabcolsep}{0pt}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{2pt}\n');
  fprintf(FID, '\\setlength{\\tabcolw}{%s-\\groupheadcolw-\\headcolw-%d\\tabcolsep}\n', par.tableWidth, 2*(nCols+1));
  fprintf(FID, '\\setlength{\\tabcolw}{\\tabcolw/%d}\n', nCols);
  fprintf(FID, '\n');
  fprintNewLength(FID, 'astwidth');
  fprintf(FID, '\\settowidth{\\astwidth}{${}^{\\ast}$}\n');
  fprintf(FID, '\\centering\n');
  fprintf(FID, '%%\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}X}\n');
  fprintf(FID, '\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}m{\\tabcolw}}\n');
  fprintf(FID, '\\newcolumntype{H}{>{\\raggedright\\arraybackslash}m{\\headcolw}}\n');
  fprintf(FID, '\\newcolumntype{G}{>{\\raggedright\\arraybackslash}m{\\groupheadcolw}}\n');
  fprintf(FID, '%%\\begin{tabularx}{\\textwidth}{ m{2cm}%s }\n', [repmat('R', 1, 2*nCols), '']); % 'rr'
  fprintf(FID, '\\begin{tabular}{ GH%s }\n', [repmat('R', 1, 2*nCols), '']); % 'rr'
  
  fprintf(FID, '\\toprule\n');
  
  % first table cells
  fprintf(FID, '{} & \\textbf{%s} & ', '$\cov$');

  % existing column groups
  if nColGroups > 0
    % header with group names
    for g = 1:nColGroups - 1
      fprintf(FID, ...
        '\\multicolumn{%d}{l}{\\parbox{%d\\tabcolw}{\\centering %s}}', ...
                   par.colGroupNum(g), par.colGroupNum(g), par.colGroups{g});
      fprintf(FID, ' & ');
    end
    fprintf(FID, ...
      '\\multicolumn{%d}{l}{\\parbox{%d\\tabcolw}{\\centering %s}}', ...
                 par.colGroupNum(nColGroups), par.colGroupNum(nColGroups), ...
                 par.colGroups{nColGroups});
    fprintf(FID, '\\\\\n');
  

    % print midrule lines
    fprintf(FID, '\\cmidrule(lr){2-2}\n');
    % group starting columns
    startColGroupVal = 3;
    colGroupStart = cumsum([startColGroupVal, par.colGroupNum]);
    for g = 1:nColGroups
      fprintf(FID, '\\cmidrule(lr){%d-%d}\n', ...
        colGroupStart(g), colGroupStart(g+1)-1 );
    end

    % header with column names
    fprintf(FID, '{%s} & ', par.rowValName);
    fprintf(FID, '{%s} & ', par.colValName);
    fprintf(FID, '{%s}', strjoin(par.colNames, '} & {')); % nCols + 1
    fprintf(FID, '\\\\\n');
    
  % non-existing column groups
  else
    % header with column names
    fprintf(FID, '{%s}', strjoin(par.colNames, '} & {')); % nCols + 1
    fprintf(FID, '\\\\\n');

    % print midrule lines
    for g = 1:nCols + 1
      fprintf(FID, '\\cmidrule(lr){%d-%d}\n', g+1, g+1);
    end

  end

  % group starting rows
  startRowGroupVal = 1;
  rowGroupStart = cumsum([startRowGroupVal, par.rowGroupNum]);
  % row loop
  for r = 1:nRows
    % row group column
    if any(r == rowGroupStart)
      if nColGroups > 0 || r > 1
        fprintf(FID, '\\midrule\n');
      end
      gId = find(r == rowGroupStart, 1, 'first');
      fprintf(FID, '\\multirow{%d}{*}{%s}', par.rowGroupNum(gId), par.rowGroups{gId});
    else
      
    end
    fprintf(FID, ' & ');
    % column with row's names
    fprintf(FID, sprintf('\\\\scriptsize{%s}', par.rowNames{r}));
    for c = 1:nCols
      % extract individual data
      dt = data(r,c);
      % mark significance
      signif = dt < par.alpha;

      fprintf(FID, ' & ');
      if signif
        fprintf(FID, '\\textbf{');
      end
      % print number
      fprintf(FID, ' \\tiny{%s}', prtTexNum(dt));
      if signif
        fprintf(FID, '}');
      end

    end
    fprintf(FID, '\\\\[-%dpt]\n', rowShrink);
  end

  % end of table
  fprintf(FID, '\\bottomrule\n');
  fprintf(FID, '%%\\end{tabularx}\n');
  fprintf(FID, '\\end{tabular}\n');

  if par.printHeader
    
    fprintf(FID, '\n');
    fprintf(FID, '\\setlength{\\tabcolsep}{\\savetabcolsep}\n');
    fprintf(FID, '\\setlength{\\cmidrulekern}{\\savecmidrulekern}\n');
    
    fprintf(FID, '\\label{tab:%s}\n', par.resultLabel);
    if par.oneColumn
      fprintf(FID, '\\end{table}\n');
    else
      fprintf(FID, '\\end{table*}\n');
    end
  end
end

function cellOfStr = formatCell(fmt, cellOfStr)
% map sprintf to each string in a cell array
  cellOfStr = cellfun(@(x) sprintf(fmt, x), cellOfStr, 'UniformOutput', false);
end

function texNum = prtTexNum(num)
% print number in exponential tex format
  NASymbol = '\makebox{---}';

  if isnan(num)
    texNum = NASymbol; %fprintf(FID, ' %s', NASymbol);
  elseif mod(num, 1) == 0
    % whole number
    texNum = num2str(num); % fprintf(FID, ' %d}', dt);
  else
    % rational number
    texNum = sprintf('%0.1e', num); %fprintf(FID, ' \\scriptsize{%0.2g}', dt);
    if numel(texNum) > 7
      texNum = sprintf('%0.0e', num); 
    end
    % remove extra zero
    if strcmp(texNum(end-2:end-1), '-0')
      texNum = texNum([1:end-2, end]);
    elseif strcmp(texNum(end-2:end-1), '+0')
      texNum = texNum([1:end-3, end]);
    end
  end
end

function fprintNewLength(FID, newLength)
% print new length to file FID
  fprintf(FID, '\\ifx\\%s\\undefined\n', newLength);
  fprintf(FID, '  \\newlength{\\%s}\n', newLength);
  fprintf(FID, '\\fi\n');
end