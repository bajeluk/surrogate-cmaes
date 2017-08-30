function [ci] = varToConfidence(yPred, sd2, alpha)
  if nargin < 3
    alpha = 0.05;
  end
  n = size(sd2, 1);
  t = tinv(1 - alpha/2, n-1);
  sd = sqrt(sd2);
  ci = [yPred - t*sd, yPred + t*sd];
end

