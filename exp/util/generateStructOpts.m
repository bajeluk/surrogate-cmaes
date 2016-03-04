function opts = generateStructOpts(carray)
% generate a struct which has all the possible value-combinations
% as defined in the field-value cellarray CARRAY
% e.g.: for carray = {'field1', { 1, 2 }, 'filed2' { 'a', 'b', 'c' } }, it will
%       produce
%     opts(1).field1 = 1; opts(1).field2 = 'a';
%     opts(2).field1 = 1; opts(2).field2 = 'b';
%     opts(3).field1 = 1; opts(3).field2 = 'c';
%     opts(4).field1 = 2; opts(4).field2 = 'a';
%     opts(5).field1 = 2; opts(5).field2 = 'b';
%     opts(6).field1 = 2; opts(6).field2 = 'c';

  lenCarray = length(carray);
  assert(mod(lenCarray, 2) == 0, 'The number of elements in field-value cell array is not even!')
  nFields    = lenCarray / 2;
  fieldNames = carray(1:2:lenCarray);
  values     = carray(2:2:lenCarray);

  % the numbers of options in each field
  % e.g.   nValues = [3 5 2]   -- |field1| = 3, |field2| = 5, |field3| = 2  (# of values)
  nValues = cellfun(@(x) length(x), values);

  % generate the final struct
  opts = struct();
  for i = 1:prod(nValues)
    multiBase = getParamIndexVector(i, nValues);
    for j = 1:nFields
      % for the j-th field, use the  'multiBase(j)'-th value from j-th values
      opts(i).(fieldNames{j}) = values{j}{multiBase(j)};
    end
  end
end
