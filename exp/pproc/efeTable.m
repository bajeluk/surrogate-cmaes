function [table, ranks] = efeTable(data, varargin)
% [table, ranks] = efeTable(data, settings)
% Creates figure of table containing rankings of expected number of 
% function evaluations (EFE) for chosen quantiles.
%
% Input:
%   data      - cell array of data
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'DataNames'   - cell array of data names (e.g. names of algorithms)
%     'DataDims'    - dimensions of data
%     'DataFuns'    - functions of data
%     'Quantiles'   - quantiles of result EFE
%     'TableDims'   - dimensions chosen to count
%     'TableFuns'   - functions chosen to count
%     'Target'      - EFE target
%
% Output:
%   table - table of EFE rankings
%   ranks - EFErankings for each function and dimension
%
% See Also:
%   rankingTable, dataReady

  % initialization
  table = [];
  if nargin < 1 || isempty(data)
    help efeTable
    return
  end
  settings = settings2struct(varargin);

  numOfData = length(data);
  datanames = defopts(settings, 'DataNames', ...
    arrayfun(@(x) ['ALG', num2str(x)], 1:numOfData, 'UniformOutput', false));
  defaultDims = [2, 3, 5, 10, 20, 40];
  funcSet.dims   = defopts(settings, 'DataDims', defaultDims(1:size(data{1}, 2)));
  funcSet.BBfunc = defopts(settings, 'DataFuns', 1:size(data{1}, 1));
  dims    = defopts(settings, 'TableDims', funcSet.dims);
  quantiles = defopts(settings, 'Quantiles', [0.25 0.5 0.75]);
  
  % create ranking table
  if ~isempty(settings)
    extraFields = {'DataNames'};
    fieldID = isfield(settings, extraFields);
    createSettings = rmfield(settings, extraFields(fieldID));
  else
    createSettings = {};
  end
  [table, ranks] = createEFETable(data, createSettings);
  
  % print table
  nEvals = length(quantiles);
  nDims = length(dims);
  maxLengthData = max([cellfun(@length, datanames), length('Quantile')]);
  colWidth = 40;
  tableSize = [11*(2+maxLengthData) + nEvals*(nDims+1)*colWidth, 20*(numOfData+2)];

  quantileRow = repmat(quantiles, [1, length(dims)+1]);
  publicTable = [quantileRow; table];
  colBase = arrayfun(@(x) repmat({[num2str(x), 'D']}, [1, nEvals]), dims, 'UniformOutput', false);
  colName = [[colBase{:}], repmat({'SUM'}, [1, nEvals])];
  rowName = [{'Quantile'}, datanames];
  f = figure('Position', [0, 0, tableSize]);
  table = uitable(f, 'Data', publicTable, ...
             'ColumnName', colName, ...
             'RowName', rowName, ...
             'ColumnWidth', {colWidth}, ...
             'Position', [1, 1, tableSize]);

end