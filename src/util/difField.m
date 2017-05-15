function [df, dVals] = difField(varargin)
% [df, dVals] = difField(s1, s2, ...) return structure different field
% values and its fieldnames.
%
% Input:
%   s1, s2, ... - structures to compare | struct or cell-array of struct
%
% Output:
%   df    - different fields among tested structures | cell-array of char
%   dVals - different values in df | cell-array
%
% Note:
%   Empty cells are ignored, empty structures are not.

  df = {};
  dVals = {};
  
  if nargin < 1
    help difField
    return
  end

  % put structrures and cell-arrays of structures to one cell-array
  structId = cellfun(@isstruct, varargin);
  cellId = cellfun(@iscell, varargin);
  assert(all(structId | cellId), 'Input is not cell-array or structure')
  
  if any(cellId)
    strCell = [varargin{cellId}];
    assert(all(cellfun(@isstruct, strCell)), 'There is a cell-array not containing a structure')
  else
    strCell = {};
  end
  strCell = [strCell, varargin{structId}];
  
  % one structure case is not comparable
  nStruct = length(strCell);
  if nStruct < 2
    return
  end
  
  % gain all subfields
  allSubfields = cellfun(@(x) subfields(x)', strCell, 'UniformOutput', false);
  uniqueSubfields = unique([allSubfields{:}]');
  nSubfields = length(uniqueSubfields);
  
  % gain field values
  subVal = cell(nSubfields, nStruct);
  for s = 1:nStruct
    for f = 1:nSubfields
      try
        subVal{f, s} = eval(['strCell{s}.', uniqueSubfields{f}]);
      catch
      end
    end
  end
  
  % gain different fields and its values
  diffID = arrayfun(@(x) ~isequal(subVal{x,:}), 1:nSubfields);
  df = uniqueSubfields(diffID);
  dVals = subVal(diffID, :);
end