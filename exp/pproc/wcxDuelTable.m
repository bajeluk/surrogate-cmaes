function wcxDuelTable(data, varargin)
% Return duelTable for two-sided Wilcoxon signed rank test's results
% in double array of size = number of dimensions x number of algorithm 
% combinations x 8 x number of FE budgets.
%
% Input:
%   data     - array of results of size = number of dimensions x number of
%              algorithm combinations x 8 x number of FE budgets | double
%            - 3rd dimension has the following 8 'columns':
%                - 1st pair member algorithm
%                - 2nd pair member algorithm
%                - number of winning instances of the 1st algorithm
%                - number of winning instances of the 2nd algorithm
%                - percentage of wins of 1st algorithm
%                - percentage of wins of 2nd algorithm
%                - p-value of Wilcoxon test
%                - p-value after Holm correction
%   settings - pairs of property (string) and value or struct with
%              properties as fields:
%
%     'Alpha'       - significance level for hypothesis testing | (0,1] |
%                     default: 0.05
%     'DataCaption' - data description inserted to table caption | string |
%                     default: 'COCO'
%     'DataDims'    - dimensions of data | integer vector | default: 2
%     'DataNames'   - cell array of data names (e.g. names of algorithms) |
%                     cell-array of string | default: {'ALG1', 'ALG2', ...}
%     'DefFile'     - name of file, where table command definitions are
%                     saved | string | default: '/tmp/wcxDefFile.tex'
%     'Evaluations' - evaluations shown in table | integer vector |
%                     default: [25, 250]
%     'HeadColW'    - header column width as a char in latex distance
%                     format | string | default: '1.7cm'
%     'HeaderFile'  - name of file with table caption and label | string |
%                     default: ''
%                   - caption with label will be written in 'ResultFile'
%                     when empty
%     'Mode'        - mode of table | {'alg', 'model'} | default: 'alg'
%     'OneColumn'   - print table to one column | boolean | default: false
%     'PrintHeader' - print table environment, caption and label | boolean
%                     | default: true
%     'ResultFile'  - file containing resulting table | string | default:
%                     '/tmp/wcxDuelTable.tex'
%     'SecondHead'  - content of heading cell of the second row (location
%                     [2, 1]) | string | default: '{\\nbFEs/\\dm}'
%     'TssList'     - list of TSS methods (mode 'model' only) | cell-array
%                     of string | default: {}
%     'TssNums'     - number of columns for each cell from TssList | double
%                     vector | default: []
%     'Vertical'    - create vertical table | boolean | default: false
%
% See Also:
%   duelTable

  
  % get basic data characteristics
  if size(data, 3) ==  8
    eightColMode = true;
    nDims = size(data, 1);
    nCombs = size(data, 2);
    nData = (1 + sqrt(1+8*nCombs))/2;
    nBudgets = size(data, 4);
  else
    eightColMode = false;
    nDims = 1;
    nData = size(data, 1);
    nBudgets = 1;
    % find p-values
    if isnumeric(varargin{1}) && all(size(data) == size(varargin{1}))
      pValData{1} = varargin{1};
      varargin = varargin(2:end);
    end
  end
  
  % parse settings
  settings = settings2struct(varargin);  
  settings.Alpha = defopts(settings, 'Alpha', 0.05);
  settings.DataCaption = defopts(settings, 'DataCaption', 'COCO');
  settings.DataDims = defopts(settings, 'DataDims', 2);
  settings.DataNames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:nData, 'UniformOutput', false));
  settings.DefFile = defopts(settings, 'DefFile', '/tmp/wcxDefFile.tex');
  settings.EightColMode = eightColMode;
  settings.Evaluations = defopts(settings, 'Evaluations', [25, 250]);
  settings.HeadColW = defopts(settings, 'HeadColW', '1.7cm');
  settings.HeaderFile = defopts(settings, 'HeaderFile', '');
  settings.Mode = defopts(settings, 'Mode', 'alg');
  settings.OneColumn = defopts(settings, 'OneColumn', false);
  settings.PrintHeader = defopts(settings, 'PrintHeader', true);
  settings.ResultFile = defopts(settings, 'ResultFile', '/tmp/wcxDuelTable.tex');
  settings.SecondHead = defopts(settings, 'SecondHead', '{\nbFEs/\dm}');
  settings.TssList = defopts(settings, 'TssList', {});
  settings.TssNums = defopts(settings, 'TssNums', []);
  settings.Vertical = defopts(settings, 'Vertical', false);
  % DEBUG mode:
  % for i = 1:numel(dataNames)
  %   fprintf('%s\n', dataNames{i})
  % end
  
  if eightColMode
    % dimension loop
    for dim = 1:nDims
      % transform results to table
      for bud = 1:nBudgets
        for com = 1:nCombs
          id1 = data(dim, com, 1, bud);
          id2 = data(dim, com, 2, bud);
          % get numbers of wins
          tableData{dim, bud}(id1, id2) = data(dim, com, 5, bud);
          tableData{dim, bud}(id2, id1) = data(dim, com, 6, bud);
          % get p-values
          pValData{dim, bud}(id1, id2) = data(dim, com, 8, bud);
          pValData{dim, bud}(id2, id1) = data(dim, com, 8, bud);
        end
      end
    end
  else
    % simple input data
    tableData{1} = data;
  end

  % print file with table definitions
  FID = fopen(settings.DefFile, 'w');
  printTexDefinitions(FID)
  fclose(FID);
  % print resulting table
  FID = fopen(settings.ResultFile, 'w');
  switch settings.Mode
    case 'alg'
      printTableTex(FID, tableData, pValData, settings)
    case 'model'
      printModelTableTex(FID, tableData, pValData, settings)
    otherwise
      error('scmaes:wcxdt:wrngmode', ...
            'Mode %s is not correct setting. Use %s instead.', ...
            settings.Mode, '''alg'' or ''model''')
  end
  fclose(FID);
end

function printTexDefinitions(FID)
% print new lengths and columns to definition file
  fprintf(FID, '\\newlength{\\dueltabcolw}\n');
  fprintf(FID, '\\newlength{\\savetabcolsep}\n');
  fprintf(FID, '\\newlength{\\savecmidrulekern}\n');
  fprintf(FID, '\\newlength{\\headcolw}\n');
  fprintf(FID, '\\newlength{\\astwidth}\n');
  
  fprintf(FID, '%%\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}X}\n');
  fprintf(FID, '\\newcolumntype{R}{>{\\raggedleft\\arraybackslash}m{\\dueltabcolw}}\n');
  fprintf(FID, '\\newcolumntype{L}{>{\\raggedright\\arraybackslash}m{\\dueltabcolw}}\n');
  fprintf(FID, '\\newcolumntype{H}{>{\\raggedright\\arraybackslash}m{\\headcolw}}\n');
end

function printTableTex(FID, table, pVals, settings)
% Prints table to file FID

  % parse settings
  alpha = settings.Alpha;
  dims = settings.DataDims;
  evaluations = settings.Evaluations;
  dataCaption = settings.DataCaption;
  datanames = settings.DataNames;
  headcolw = settings.HeadColW;
  headerFile = settings.HeaderFile;
  oneColumn = settings.OneColumn;
  printHeader = settings.PrintHeader;
  secondHead = settings.SecondHead;
  vertical = settings.Vertical;
  
  numOfData = length(datanames);

  % in case of multiple dimensions in one table insert char
  if iscell(dims)
    dimCell = dims;
    dims = 0;
  else
    dimCell = '';
  end
  nEvals = length(evaluations);
  nDims = length(dims);
  
  % symbol for number of evaluations reaching the best target
  bestSymbol = '\bestFED';
  maxFunEvalsSymbol = '\maxFED';
  maxFunEvalsString = '$250\dm$';
  ftargetString = '10^{-8}';
  NASymbol = '\makebox{---}';

  % representation of evaluation counts as a fraction
  evaluationsString = cell(1, nEvals);
  for e = 1:nEvals
    s = strsplit(strtrim(rats(evaluations(e))), '/');
    if length(s) == 1
      evaluationsString{e} = sprintf('%s\\\\mbox{\\\\hspace{\\\\astwidth}}', s{1});
    else
      evaluationsString{e} = sprintf('{\\\\large\\\\sfrac{%s}{%s}}\\\\mbox{\\\\hspace{\\\\astwidth}}', s{1}, s{2});
    end
  end
  
  % header init
  if printHeader
    % checkout where to print header
    if isempty(headerFile)
      headerFID = FID;
      if vertical
        fprintf(FID, '\\begin{sidewaystable}\n');
      elseif oneColumn
        fprintf(FID, '\\begin{table}[t]\n');
      else
        fprintf(FID, '\\begin{table*}[t]\n');
      end
    % header in extra file
    else
      headerFID = fopen(headerFile, 'w');
    end
  end

  fprintf(FID, '\\setlength{\\savetabcolsep}{\\tabcolsep}\n');
  fprintf(FID, '\\setlength{\\savecmidrulekern}{\\cmidrulekern}\n');
  fprintf(FID, '\n');
  fprintf(FID, '\\setlength{\\headcolw}{%s}\n', headcolw);
  fprintf(FID, '\\setlength{\\tabcolsep}{0pt}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{2pt}\n');
  fprintf(FID, '\\setlength{\\dueltabcolw}{%0.1f\\textwidth-\\headcolw-%d\\tabcolsep}\n', ...
               1-0.6*oneColumn, nEvals*(nEvals*numOfData+1));
  fprintf(FID, '\\setlength{\\dueltabcolw}{\\dueltabcolw/%d}\n', nEvals*numOfData);
  fprintf(FID, '\n');
  
  fprintf(FID, '\\settowidth{\\astwidth}{${}^{\\ast}$}\n');
  fprintf(FID, '\\centering\n');
  
  if printHeader
    % caption printing
    if isempty(dimCell)
      dimString = sprintf('$%d\\\\dm$', dims(1));
      for dim = 2:nDims-1
        dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dims(dim))}, ', ');
      end

      if nDims > 1
        dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dims(nDims))}, ' and ');
      end
    else
      dimString = sprintf('$%d\\\\dm$', dimCell{1});
      for dim = 2:numel(dimCell)-1
        dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dimCell{dim})}, ', ');
      end

      if numel(dimCell) > 1
        dimString = strjoin({dimString, sprintf('$%d\\\\dm$', dimCell{end})}, ' and ');
      end
    end
    
    fprintf(headerFID, ...
      ['\\caption{A pairwise comparison of the evolution controls, ', ...
       'models, and their combinations in ', dimString, ...
       ' over the %s for different evaluation budgets.\n', ...
       'The percentage of wins of $i$-th algorithm against $j$-th algorithm ', ...
       'over all benchmark instances is given in the $i$-th row and $j$-th column.\n'], ...
      dataCaption);
    
    fprintf(headerFID, ['The numbers in bold mark the row algorithm being ', ...
      'significantly better than the column algorithm ', ...
      'according to the two-sided Wilcoxon signed rank test with the Holm ', ...
      'correction at family-wise significance level $\\alpha=%.2f$.\n'], ...
      alpha);
    
    fprintf(headerFID, ['%%%s\\ denotes the smallest \\nbFEs\\ at which any of the tested algorithms ', ...
      'reached the target $\\ftarget = %s$.\n', ...
      '%%%s\\ is equal to %s, the overall budget for the experiments.\n', ...
      '}\n'], ...
      bestSymbol, ftargetString, maxFunEvalsSymbol, maxFunEvalsString ...
      );
    
    fprintf(headerFID, '\n');
    [~, labelTab] = fileparts(settings.ResultFile);
    fprintf(headerFID, '\\label{tab:%s}\n', labelTab);
  end
  
  % uncomment the following line and comment the next for tabularx package
  % usage
  % fprintf(FID, '\\begin{tabularx}{\\textwidth}{ m{2cm}%s }\n', [repmat('R', 1, nEvals*numOfData), '']); % 'rr'
  fprintf(FID, '\\begin{tabular}{ H%s }\n', [repmat('L', 1, nEvals*numOfData), '']); % 'rr'
  
  for dim = 1:nDims
    if dim == 1
      fprintf(FID, '\\toprule\n');
    else
      fprintf(FID, '\\midrule[\\heavyrulewidth]\n');
    end

    if ~isempty(dimCell)
      fprintf(FID, '$\\bm{%d-%d\\dm}$ &', dimCell{1}, dimCell{end});
    else
      fprintf(FID, '$\\bm{%d\\dm}$ &', dims(dim));
    end

    % header with algorithm names
    mcLine = ['\\\\multicolumn{', num2str(nEvals), '}{l}{\\\\parbox{', num2str(nEvals), '\\\\dueltabcolw}{\\\\centering %s}}'];
    fprintf(FID, strjoin(formatCell(mcLine, datanames), ' & '));
    fprintf(FID, '\\\\\n');

    fprintf(FID, '\\cmidrule(lr){1-1}\n');
    for i = 1:numOfData
      fprintf(FID, '\\cmidrule(lr){%d-%d}\n', nEvals*(i-1)+2, nEvals*i+1);
    end

    % header with evaluation numbers
    fprintf(FID, '%s & ', secondHead);
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
            % mark significance according to a powerful correction by one
            % star
            sigdiff = dt(i, j) > dt(j, i) && pv(i, j) < alpha;

            fprintf(FID, ' & ');
            % show significance in bold
            if sigdiff
              fprintf(FID, '\\textbf{');
            end
            % print the number
            if mod(dt(i, j), 1) == 0
              % whole number
              fprintf(FID, '%d', dt(i, j));
            else
              % rational number
              fprintf(FID, '%0.1f', dt(i, j));
            end
            % finish significance
            if sigdiff
              fprintf(FID, '}');
            end
          else
            % symbols at diagonal
            fprintf(FID, ' & %s', NASymbol);
          end
        end
      end
      fprintf(FID, '\\\\\n');
    end

    if dim == nDims
      fprintf(FID, '\\bottomrule\n');
    end
  end
  % uncomment the following line and comment the next for tabularx package
  % usage
  % fprintf(FID, '\\end{tabularx}\n');
  fprintf(FID, '\\end{tabular}\n');

  fprintf(FID, '\\setlength{\\tabcolsep}{\\savetabcolsep}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{\\savecmidrulekern}\n');
  if printHeader
    if isempty(headerFile)
      if vertical
        fprintf(FID, '\\end{sidewaystable}\n');
      elseif oneColumn
        fprintf(FID, '\\end{table}\n');
      else
        fprintf(FID, '\\end{table*}\n');
      end
    else
    % close header file
      fclose(headerFID);
    end
  end
end

function printModelTableTex(FID, table, pVals, settings)
% Prints table to file FID

  % parse settings
  alpha = settings.Alpha;
  models = settings.Evaluations;
  dataCaption = settings.DataCaption;
  datanames = settings.DataNames;
  headcolw = settings.HeadColW;
  headerFile = settings.HeaderFile;
  oneColumn = settings.OneColumn;
  printHeader = settings.PrintHeader;
  tssList = settings.TssList;
  tssNums = settings.TssNums;
  vertical = settings.Vertical;

  % definitions
  nModelSettings = [8, 9, 1, 1, 8, 9, 1, 1, 1];
  nModels = length(models);

  % symbol definitions
  NASymbol = ''; % '\makebox{---}';

  % header init
  if printHeader
    % checkout where to print header
    if isempty(headerFile)
      headerFID = FID;
      if vertical
        fprintf(FID, '\\begin{sidewaystable}\n');
      elseif oneColumn
        fprintf(FID, '\\begin{table}[t]\n');
      else
        fprintf(FID, '\\begin{table*}[t]\n');
      end
    % header in extra file
    else
      headerFID = fopen(headerFile, 'w');
    end
  end

  fprintf(FID, '\\setlength{\\savetabcolsep}{\\tabcolsep}\n');
  fprintf(FID, '\\setlength{\\savecmidrulekern}{\\cmidrulekern}\n');
  fprintf(FID, '\n');
  fprintf(FID, '\\setlength{\\headcolw}{%s}\n', headcolw);
  fprintf(FID, '\\setlength{\\tabcolsep}{0pt}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{2pt}\n');
  fprintf(FID, '\\setlength{\\dueltabcolw}{%0.1f\\textwidth-\\headcolw-%d\\tabcolsep}\n', ...
               1.3-0.6*oneColumn, nModels*(nModels+1));
  fprintf(FID, '\\setlength{\\dueltabcolw}{\\dueltabcolw/%d}\n', nModels);
  fprintf(FID, '\n');

  fprintf(FID, '\\settowidth{\\astwidth}{${}^{\\ast}$}\n');
  fprintf(FID, '\\centering\n');

  if printHeader
    % caption printing
    fprintf(headerFID, '\\caption{\n');
    fprintf(headerFID, ['  A pairwise comparison of the model settings %s in different TSS.\n', ...
       'The percentage of wins of $i$-th model setting against $j$-th model setting ', ...
       'over all available data is given in the $i$-th row and $j$-th column.\n'], ...
      dataCaption);

    fprintf(headerFID, ['  The numbers in bold mark the row model setting being ', ...
      'significantly better than the column model setting ', ...
      'according to the two-sided Wilcoxon signed rank test with the Holm ', ...
      'correction at family-wise significance level $\\alpha=%.2f$.\n'], ...
      alpha);
    fprintf(headerFID, '}\n');
    fprintf(headerFID, '\n');
    [~, labelTab] = fileparts(settings.ResultFile);
    fprintf(headerFID, '\\label{tab:%s}\n', labelTab);
  end

  % resize to textwidth if necessary
  fprintf(FID, '\\resizebox{\\textwidth}{!}{%%\n');
  % uncomment the following line and comment the next for tabularx package
  % usage
  % fprintf(FID, '\\begin{tabularx}{\\textwidth}{ m{2cm}%s }\n', [repmat('R', 1, nModels), '']); % 'rr'
  fprintf(FID, '\\begin{tabular}{ HHH%s }\n', [repmat('L', 1, nModels), '']); % 'rr'
  % first top (thick) line
  fprintf(FID, '\\toprule\n');
  % first three cells
  fprintf(FID, '\\multicolumn{2}{c}{\\multirow{2}{*}{\\textbf{%s}}} & TSS & ', dataCaption);

  % header with TSS
  tssLine = @(x, y) sprintf('\\\\multicolumn{%d}{l}{\\\\parbox{%d\\\\dueltabcolw}{\\\\centering %s}}', x, x, y{1});
  fprintf(FID, strjoin(arrayfun(tssLine, tssNums, tssList, 'Uni', false), ' & '));
  fprintf(FID, '\\\\\n');

  % data columns mid-lines
  cmidruleCols = [0, cumsum(tssNums), sum(tssNums)];
  for i = 2:numel(tssNums)+1
    fprintf(FID, '\\cmidrule(lr){%d-%d}\n', cmidruleCols(i-1)+4, cmidruleCols(i)+3);
  end

  % header with model names
  fprintf(FID, '\\multicolumn{2}{c}{} & model & ');
  mcLine = @(x, y) sprintf('\\\\multicolumn{%d}{l}{\\\\parbox{%d\\\\dueltabcolw}{\\\\centering %s}}', x, x, y{1});
  fprintf(FID, strjoin(arrayfun(mcLine, nModelSettings, datanames, 'Uni', false), ' & '));
  fprintf(FID, '\\\\\n');

  % data columns mid-lines
  cmidruleCols = [0, cumsum(nModelSettings), sum(nModelSettings)];
  for i = 2:numel(nModelSettings)+1
    fprintf(FID, '\\cmidrule(lr){%d-%d}\n', cmidruleCols(i-1)+4, cmidruleCols(i)+3);
  end

  % header with model settings
  fprintf(FID, 'TSS & model & setting & ');
  fprintf(FID, strjoin(models, ' & ')); % numOfData + 1
  fprintf(FID, '\\\\\n');
  fprintf(FID, '\\midrule\n');

  for i = 1:sum(nModelSettings)
    % first column (TSS)
    tssId = i == [1, cumsum(tssNums(1:end-1)) + 1];
    if any(tssId)
      fprintf(FID, sprintf('\\\\multirow{%d}{*}{%s}', ...
                           tssNums(tssId), tssList{tssId}));
    end
    fprintf(FID, ' & ');

    % second column (model type)
    mdlId = i == [1, cumsum(nModelSettings(1:end-1)) + 1];
    if any(mdlId)
      fprintf(FID, sprintf('\\\\multirow{%d}{*}{%s}', ...
                           nModelSettings(mdlId), datanames{mdlId}));
    end
    fprintf(FID, ' & ');

    % column with model settings's name
    fprintf(FID, sprintf('%s', models{i}));

    for j = 1:sum(nModelSettings)
      dt = table{1};
      pv = pVals{1};

      if i ~= j
        % mark significance according to a powerful correction by one star
        sigdiff = dt(i, j) > dt(j, i) && pv(i, j) < alpha;

        fprintf(FID, ' & ');
        % print color background
        fprintf(FID, '\\cellcolor[rgb]{1, %0.3f, %0.3f} ', ...
                1 - dt(i, j)/100, (1 - dt(i,j)/100)*0.6 + 0.4);
        % show significance in bold
        if sigdiff
          fprintf(FID, '\\textbf{');
        end
        % print the number
        if mod(dt(i, j), 1) == 0
          % whole number
          fprintf(FID, '%d', dt(i, j));
        else
          % rational number
          fprintf(FID, '%0.1f', dt(i, j));
        end
        % finish significance
        if sigdiff
          fprintf(FID, '}');
        end
      else
        % symbols at diagonal
        fprintf(FID, ' & %s', NASymbol);
      end
    end
    % end of line
    fprintf(FID, '\\\\\n');

    % split TSS lines
    if any(i == cumsum(tssNums(1:end-1)))
      fprintf(FID, '\\cmidrule(lr){1-3}\n');
    % split model type lines
    % elseif any(i == cumsum(nModelSettings(1:end-1)))
    %   fprintf(FID, '\\cmidrule(lr){2-3}\n');
    end
  end

  % last line
  fprintf(FID, '\\bottomrule\n');
  % end of tabular environment
  % uncomment the following line and comment the next for tabularx package
  % usage
  % fprintf(FID, '\\end{tabularx}\n');
  fprintf(FID, '\\end{tabular}\n');
  % end resizebox environment
  fprintf(FID, '}\n');

  fprintf(FID, '\\setlength{\\tabcolsep}{\\savetabcolsep}\n');
  fprintf(FID, '\\setlength{\\cmidrulekern}{\\savecmidrulekern}\n');
  if printHeader
    if isempty(headerFile)
      if vertical
        fprintf(FID, '\\end{sidewaystable}\n');
      elseif oneColumn
        fprintf(FID, '\\end{table}\n');
      else
        fprintf(FID, '\\end{table*}\n');
      end
    else
    % close header file
      fclose(headerFID);
    end
  end
end

function cellOfStr = formatCell(fmt, cellOfStr)
% map sprintf to each string in a cell array
  cellOfStr = cellfun(@(x) sprintf(fmt, x), cellOfStr, 'UniformOutput', false);
end