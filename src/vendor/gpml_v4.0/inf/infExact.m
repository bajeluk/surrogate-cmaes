function [post nlZ dnlZ] = infExact(hyp, mean, cov, lik, x, y)
  [post nlZ dnlZ] = infGaussLik(hyp, mean, cov, lik, x, y);
end
