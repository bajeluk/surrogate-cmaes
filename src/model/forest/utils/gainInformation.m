function [gain] = gainInformation(current, left, right)
  d = size(current.y, 2);
  if isfield(current.sigma)
    % covariance matrix
    entropyFunc = @(y)
    current.entropy = 
  else
    % only variance
  end
  current.entropy = 
  entropyCurrent = current.sd2;
  entropyLeft = left.sd2
  entropyRight = immse(right.y, right.yPred);
  entropyNew = (entropyLeft * numel(left.y) + entropyRight * numel(right.y)) ...
    / (numel(left.y) + numel(right.y));
  gain = entropyCurrent - entropyNew;
end