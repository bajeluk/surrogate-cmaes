function sc = combineFieldValues(s)
% sc = combineFieldValues(s) creates cell-array of structures containing 
% all combinations of the input structure (or cell-array of structures) 
% with cell-arrays as fields.
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
  
  if ~iscell(s)
    s = {s};
  end
  
  sc = {};
  nsc = 0;
  for j = 1 : numel(s)
    sFields = fieldnames(s{j});
    nFields = length(sFields);

    % check input has cell-array fields
    isNotCellField = find(~cellfun(@(x) iscell(s{j}.(x)), sFields));
    % not cell-array fields are changed to cells
    for f = 1:length(isNotCellField)
      s{j}.(sFields{isNotCellField(f)}) = {s{j}.(sFields{isNotCellField(f)})};
    end

    % number of settings in fields
    nVals = cellfun(@(x) length(s{j}.(x)), sFields);
    sc(end+1:end+prod(nVals)) = cell(prod(nVals), 1);
    for i = 1:prod(nVals)
      parId = getParamIndexVector(i, nVals);
      nsc = nsc + 1;
      for f = 1:nFields
        sc{nsc}.(sFields{f}) = s{j}.(sFields{f}){parId(f)};
      end
    end
  end
  
end