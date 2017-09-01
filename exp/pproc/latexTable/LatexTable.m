classdef LatexTable < handle
% An easy to use class that generates a LaTeX table from a given MATLAB
% cellarray or table.
%
% Author:       Lukas Bajer
% Original author:  Eli Duenisch
% Date:         July 31, 2017
% License:      This code is licensed using BSD 2 to maximize your freedom of using it :)
%
% original copyright:
% ----------------------------------------------------------------------------------
%  Copyright (c) 2016, Eli Duenisch
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%  
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
%  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ----------------------------------------------------------------------------------
%
% % Formatting-string to set the precision of the table values:
% % For using different formats in different rows use a cell array like
% % {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% % where myFormatString_ are formatting-strings and numberOfValues_ are the
% % number of table columns or rows that the preceding formatting-string applies.
% % Please make sure the sum of numberOfValues_ matches the number of columns or
% % rows in input.tableData!
% %
% % opts.dataFormat = {'%.4f'}; % uses three digit precision floating point for all data values
% opts.dataFormat = {'%.3f',2,'%.1f',1}; % three digits precision for first two columns, one digit for the last
%
% % Define how NaN values in opts.tableData should be printed in the LaTex table:
% opts.dataNanString = '$-$';
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% opts.tableColumnAlignment = 'c';
% % OR
% opts.tableColumnAlignment = {'c', 'l', 'r', 'c'};
%
% % Switch table borders on/off:
% opts.tableBorders = 1;
%
% % Switch table booktabs on/off:
% opts.booktabs = 1;
%
% % LaTex table caption:
% opts.tableCaption = 'MyTableCaption';
%
% % LaTex table label:
% opts.tableLabel = 'MyTableLabel';
%
%
% header = {'covf', 'meanf', 'ell', 'D2', 'D3', 'D5', 'D10', 'D20', 'average'};
% data   = [covCol, meanCol, ellCol, num2cell(rdePerDim), meanRDECol];
% 
% fixedHypersTable = cell2table(data, 'VariableNames', header);
% writetable(fixedHypersTable, '../latex_scmaes/ec2016paper/data/fixedHypers.csv');
% 
% lt = LatexTable(fixedHypersTable);
% lt.headerRow = {'covariance f.', '$m_\mu$', '$\ell$', '$2D$', '$3D$', '$5D$', '$10D$', '$20D$', '\textbf{average}'}';
% lt.opts.tableColumnAlignment = num2cell('lcccccccc');
% lt.opts.numericFormat = '$%.2f$';  % set up the default format for numerics
% lt.opts.booktabs = 1;     % use nice-looking tables
% lt.opts.latexHeader = 0;  % do not print the header like '\begin{table}...
% % identify minimas and set them bold
% [~, minRows] = min([rdePerDim, meanRDE]);
% for j = 1:size([rdePerDim, meanRDE], 2)
%   lt.setFormatXY(minRows(j), 3+j, '$\\bf %.2f$');
% end
% latexRows = lt.toStringRows(lt.toStringTable);
% % delete the lines \begin{tabular}{...} \toprule
% % and              \bottomrule  \end{tabular}
% latexRows([1,2,end-1,end]) = [];
% % save the result in the file
% fid = fopen('../latex_scmaes/ec2016paper/data/fixedHypers.tex', 'w');
% for i = 1:length(latexRows)
%   fprintf(fid, '%s\n', latexRows{i});
% end
% fclose(fid);

  properties
    opts
    data
    headerRow
    headerCol
    formats
  end
  
  methods
    function this = LatexTable(data, opts)
      %%%%%%%%%%%%%%%%%%%%%%%%%% Default settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % These settings are used if the corresponding optional inputs are not given.
      %
      % Placement of the table in LaTex document
      if (~exist('opts', 'var') || ~isstruct(opts))
        fprintf('Options are not specified...');
        opts = struct();
      end

      if (isfield(opts,'tablePlacement') && ~isempty(opts.tablePlacement))
        opts.tablePlacement = ['[', opts.tablePlacement, ']'];
      else
        opts.tablePlacement = '';
      end

      this.opts = opts;

      % Pivoting of the input data switched off per default:
      % Sets the default display format of numeric values in the LaTeX table to '%.4f'
      % (4 digits floating point precision).
      if ~isfield(opts,'numericFormat'), this.opts.numericFormat = '%.4f'; end
      % Define what should happen with NaN values in opts.tableData:
      if ~isfield(opts,'dataNanString'), this.opts.dataNanString = '$-$'; end
      % Specify the alignment of the columns:
      % 'l' for left-justified, 'c' for centered, 'r' for right-justified
      if ~isfield(opts,'tableColumnAlignment'), this.opts.tableColumnAlignment = 'c'; end
      % Print also LaTeX header like '\begin{table}...\caption{...}..\begin{table}
      if ~isfield(opts,'latexHeader'), this.opts.latexHeader = 1; end
      % Specify whether the table has borders:
      % 0 for no borders, 1 for borders
      if ~isfield(opts,'tableBorders'), this.opts.tableBorders = 0; end
      % Specify whether to use booktabs formatting or regular table formatting:
      if ~isfield(opts,'booktabs')
        this.opts.booktabs = 0;
      else
        if (opts.booktabs)
          this.opts.tableBorders = 0;
        end
      end
      % Other optional fields:
      if ~isfield(opts,'tableCaption'), this.opts.tableCaption = ''; end
      if ~isfield(opts,'tableLabel'), this.opts.tableLabel = ''; end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      this.headerRow = {};
      this.headerCol = {};

      % convert the input 'data' into cell array
      if isa(data, 'table')
        if(~isempty(data.Properties.RowNames))
          this.headerCol = data.Properties.RowNames';
        end
        if(~isempty(data.Properties.VariableNames))
          this.headerRow = data.Properties.VariableNames';
        end
        this.data = table2cell(data);
      elseif (isnumeric(data))
        this.data = num2cell(data);
      elseif (iscell(data))
        this.data = data;
      else
        error('Data has to be either table, cell array, or numeric matrix');
      end
      this.formats = cell(size(this.data));
    end
    
    function n = getNRows(this)
      % get number of rows
      n = size(this.data, 1);
    end
    
    function n = getNCols(this)
      % get number of columns
      n = size(this.data, 2);
    end

    function setHeaderRow(this, header)
      % set the following nx1 cell array of strings as header-row of the table
      this.headerRow = header;
    end
    
    function setHeaderCol(this, colHeader)
      % set the following rx1 cell array of strings as header-row of the table
      this.headerCol = colHeader;
    end

    function setColumnFormat(this, colFormats)
      % specify the 'printf'-like format for the columns, or for all
      % columns
      if (ischar(colFormats))
        colFormats = {colFormats};
      elseif (~iscell(colFormats))
        error('Formats has to be in a cell array');
      end
      
      if (length(colFormats) == 1)
        colFormats = repmat(colFormats, this.getNRows(), this.getNCols());
      elseif (length(colFormats) ~= this.getNCols())
        error('The number of cells in colFormats does not agree with # cols of the data');
      end
      
      this.formats = repmat(colFormats, this.nRows, 1);
    end
    
    function setFormatXY(this, i, j, formatStr)
      % set the printf-like format for the coordinates x,y
      if (~ischar(formatStr))
        error('Format has to be a string');
      end
      this.formats{i,j} = formatStr;
    end
    
    function prependFormatXY(this, i, j, formatStr)
      % prepend the printf-like format at the coordinates x,y
      if (isempty(this.formats{i,j}))
        if (isnumeric(this.data{i,j}))
          this.formats{i,j} = this.opts.numericFormat;
        else
          this.formats{i,j} = '%s';
        end
      end
      this.formats{i,j} = [formatStr, ' ', this.formats{i,j}];
    end

    function appendFormatXY(this, i, j, formatStr)
      % prepend the printf-like format at the coordinates x,y
      if (isempty(this.formats{i,j}))
        if (isnumeric(this.data{i,j}))
          this.formats{i,j} = this.opts.numericFormat;
        else
          this.formats{i,j} = '%s';
        end
      end
      this.formats{i,j} = [this.formats{i,j}, ' ', formatStr];
    end
    
    function colorizeSubMatrixInGray(this, values, row, col, minGray, maxGray)
      minValue = min( min( values ) );
      maxValue = max( max( values ) );
      for i = 1:size(values,1)
        for j = 1:size(values,2)
          grayValue = ((values(i,j) - minValue) / (maxValue-minValue)) ...
              * (maxGray - minGray) + minGray;
          if (grayValue < 1.0)
            this.prependFormatXY(row+i-1, col+j-1, ['\\cellcolor[gray]{' sprintf('%.4f', grayValue) '}']);
          end
        end
      end
    end

    function stringTable = toStringTable(this)
      % generate a cell array of strings wich would go into the final table
      % (after separating them with '&')
      stringTable = {};
      isHeaderRow = 0;
      isHeaderCol = 0;

      if (~isempty(this.headerCol))
        isHeaderCol = 1;
      end
      if (~isempty(this.headerRow))
        stringTable = this.headerRow;
        isHeaderRow = 1;
      end

      for row = 1:(this.getNRows())
        if (isHeaderCol)
          stringTable{row+isHeaderRow, 1} = this.headerCol{row};
        end
        
        for col = 1:(this.getNCols())
          dataValue = this.data{row, col};
          if isnan(dataValue)
            dataValue = this.opts.dataNanString;
          elseif (isnumeric(dataValue) || ischar(dataValue))
            if (~isempty(this.formats{row, col}))
              thisFormat = this.formats{row, col};
            elseif (isnumeric(dataValue))
              thisFormat = this.opts.numericFormat;
            else
              thisFormat = '%s';
            end
            if (~isempty(dataValue))
              dataValue = sprintf(thisFormat, dataValue);
            else
              dataValue = '';
            end
          end
          stringTable{row+isHeaderRow, isHeaderCol + col} = dataValue;
        end
      end
    end
    
    function latex = toStringRows(this, stringTable)
      % from the cell array of final LaTeX strings, generate the
      % final LaTeX code into 'latex' -- cell array of lines of code
      %
      % make table header lines:
      hLine = '\hline';
      latex = {};
      if (this.opts.latexHeader)
        latex = {['\begin{table}',this.opts.tablePlacement]; '\centering'};
        if (~isempty(this.opts.tableCaption))
          latex{end+1} = ['\caption{',this.opts.tableCaption,'}'];
        end
        if (~isempty(this.opts.tableLabel))
          latex{end+1} = ['\label{table:',this.opts.tableLabel,'}'];
        end
      end
      % set up the alignment string
      if (ischar(this.opts.tableColumnAlignment))
        this.opts.tableColumnAlignment = {this.opts.tableColumnAlignment};
      end
      if (length(this.opts.tableColumnAlignment) == 1)
        align = repmat(this.opts.tableColumnAlignment, 1, this.getNCols());
      else
        align = this.opts.tableColumnAlignment;
      end
      if (this.opts.tableBorders)
        header = ['\begin{tabular}','{|',strjoin(align, '|'),'|}'];
      else
        header = ['\begin{tabular}','{',strjoin(align, ''),'}'];
      end
      latex{end+1} = header;

      % generate table itself
      if (this.opts.booktabs)
        latex(end+1) = {'\toprule'};
      else
        if (this.opts.tableBorders)
          latex(end+1) = {hLine};
        end
      end

      % identify the maximal widths of columns for the right text alignment
      % of the '&' characters
      colWidths = zeros(size(stringTable, 2));
      for j = 1:size(stringTable, 2)
        colWidths(j) = 0;
        for i = 1:size(stringTable, 1)
          colWidths(j) = max(colWidths(j), length(stringTable{i,j}));
        end
      end

      % generate the table line-by-line
      for i = 1:size(stringTable, 1)
        thisRow = '';
        skip = 0;
        nCols = size(stringTable, 2);
        
        for j = 1:nCols
          [~, multiCols] = regexp(stringTable{i,j}, '\\multicolumn{([0-9]+)}', 'match', 'tokens');
          thisRow = [thisRow sprintf(['%' num2str(colWidths(j)) 's'], stringTable{i,j})];
          if (~isempty(multiCols) && str2num(multiCols{1}{1}) > 1)
            skip = skip + (str2num(multiCols{1}{1}) - 1);
            thisRow = [thisRow, '   '];
          elseif (j < nCols)
            thisRow = [thisRow, ' & '];
          end
          % if (j >= nCols - skip)
          %   break;
          % end
        end
        latex{end+1} = [thisRow ' \\'];

        if (i == 1)
          if (this.opts.booktabs)
            latex{end+1} = '\midrule';
          elseif (this.opts.tableBorders)
            latex(end+1) = {hLine};
          end
        end
      end
 
      if (this.opts.booktabs)
        latex(end+1) = {'\bottomrule'};
      elseif (this.opts.tableBorders)
        latex(end+1) = {hLine};
      end
      latex(end+1) = {'\end{tabular}'};

      if (this.opts.latexHeader)
        % make footer lines for table:
        latex(end+1) = {'\end{table}'};
      end
    end
  end
end
