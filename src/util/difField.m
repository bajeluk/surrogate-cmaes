function [df, dVals] = difField(struct1, struct2)
% return structure different field values and its fieldnames

  % put structrures and cell-arrays of structures to one cell-array
  if isstruct(struct1) 
    strCell{1} = struct1;
  elseif iscell(struct1)
    strCell = struct1;
  else
    error('struct1 is not cell-array or structure')
  end
  if nargin == 2
    if isstruct(struct2)
      strCell{end+1} = struct2; 
    elseif iscell(struct2)
      strCell(end+1:end+length(struct2))
    else
      error('struct2 is not cell-array or structure')
    end
  end
  
  df = {};
  dVals = {};
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