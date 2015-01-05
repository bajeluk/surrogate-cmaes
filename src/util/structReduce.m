function out = structReduce(structin, func, varargin)
% map primitive implementation for struct array

if nargin > 2
  out = varargin{1};
else
  out = [];
end

for i = 1:length(structin)
  out = feval(func, out, structin(i));
end
