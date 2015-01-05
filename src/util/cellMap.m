function cout = cellMap(cin, func)
% map primitive implementation for cell array
cout = cell(size(cin));

for i = 1:length(cin)
  cout{i} = feval(func, cin{i});
end
