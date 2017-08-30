function [sd2] = confidenceToVar(ci, alpha)
  if nargin < 2
    alpha = 0.05;
  end
  n = size(ci, 1);
  t = tinv(1 - alpha/2, n-1);
  sd = (ci(:, 2) - ci(:, 1)) / (2*t);
  sd2 = sd.^2;
end

