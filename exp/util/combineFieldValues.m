function sc = combineFieldValues(s)
% sc = combineFieldValues(s) creates cell-array of structures containing 
% all combinations of the input structure with cell-arrays as fields.
%
% Example:
%   s = struct('a', {{1,2,3}}, 'b', {{'c','d'}})
%   s = 
%        a: {[1]  [2]  [3]}
%        b: {'c'  'd'}
%
%   sc = combineFieldValues(s);
%   sc{1} = 
%            a: 1
%            b: 'c'
%   sc{6} = 
%            a: 3
%            b: 'd'
%
% See Also:
%   getParamIndexVector
  
  sFields = fieldnames(s);
  nFields = length(sFields);
  
  % check input has cell-array fields
  isNotCellField = find(~cellfun(@(x) iscell(s.(x)), sFields));
  % not cell-array fields are changed to cells
  for f = 1:length(isNotCellField)
    s.(sFields{isNotCellField(f)}) = {s.(sFields{isNotCellField(f)})};
  end
  
  % number of settings in fields
  nVals = cellfun(@(x) length(s.(x)), sFields);
  sc = cell(prod(nVals), 1);
  for i = 1:prod(nVals)
    parId = getParamIndexVector(i, nVals);
    for f = 1:nFields
      sc{i}.(sFields{f}) = s.(sFields{f}){parId(f)};
    end
  end
  
end