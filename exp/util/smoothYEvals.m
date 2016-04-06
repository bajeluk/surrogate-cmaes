function [yopt, evals] = smoothYEvals(y_evals, N)

c = 1;
y_evals = [NaN, 1; y_evals];
curr_min = NaN;

if (nargin < 2)
  N = 500;
end

evals = (1:N)';
yopt  = zeros(N,1);
for i = 1:N
  while (size(y_evals,1) > c  &&  y_evals(c+1,2) == i)
    c = c + 1;
    curr_min = min(curr_min, y_evals(c,1));
  end
  yopt(i) = curr_min;
end
