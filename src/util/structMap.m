function out = structMap(structin, func)
% map primitive implementation for struct array
out = cell(size(structin));

for i = 1:length(structin)
  % feval(func, structin(i));
  out{i} = feval(func, structin(i));
end
