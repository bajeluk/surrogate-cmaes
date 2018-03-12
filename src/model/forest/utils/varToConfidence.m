function ci = varToConfidence(yPred, sd2, alpha)
% ci = varToConfidence(yPred, sd2, alpha) returns confidence interval for
% given predicted values yPred and variance sd2 on level alpha.

  if nargin < 3
    if nargin < 1
      help varToConfidence
      if nargout > 0
        ci = [];
      end
      return
    end
    alpha = 0.05;
  end
  
  n = size(sd2, 1);
  t = tinv(1 - alpha/2, n-1);
  sd = sqrt(sd2);
  ci = [yPred - t*sd, yPred + t*sd];
end

