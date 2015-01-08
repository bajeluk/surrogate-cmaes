function [yopt evals] = smoothStatEvals(y_evals, varargin)

c = 1;
y_evals = [NaN 1; y_evals];
curr_value = NaN;

if (nargin > 1)
  N = varargin{1};
else
  N = 500;
end

evals = 1:N';
yopt  = zeros(N,1);
for i = 1:N
  while (size(y_evals,1) > c  &&  y_evals(c+1,2) == i)
    c = c + 1;
    curr_value = y_evals(c,1);
  end
  yopt(i) = curr_value;
end
