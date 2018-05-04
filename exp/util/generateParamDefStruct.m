function defs = generateParamDefStruct(carray)
  assert(mod(length(carray), 2) == 0, 'The number of elements in field-value cell array is not even!')
  defs = struct();
  for i = 1:(length(carray)/2)
    defs(i).name   = carray{(i-1)*2 + 1};
    defs(i).values = carray{i*2};
  end
end