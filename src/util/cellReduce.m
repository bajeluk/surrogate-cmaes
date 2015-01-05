function out = cellReduce(cin, func, varargin)
% map primitive implementation for cell array

if nargin > 2
  out = varargin{1};
else
  out = [];
end

for i = 1:length(cin)
  out = feval(func, out, cin{i});
end
