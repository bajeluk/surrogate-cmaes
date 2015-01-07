function [medians q1 q3 evals] = statisticsFromYEvals(y_evals, varargin)
  if (nargin > 1)
    N = varargin{1};
  else
    N = 500;
  end

  evals = 1:N';

  M = zeros(N,length(y_evals));
  for i = 1:length(y_evals)
    M(:,i) = smoothYEvals(y_evals{i}, N);
  end
  Q = quantile(M, [0.25 0.5 0.75], 2);
  q1 = Q(:,1);
  medians = Q(:,2);
  q3 = Q(:,3);
end
