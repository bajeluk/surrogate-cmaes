function [gain] = gainMse(current, left, right)
  current.mse = immse(current.y, current.yPred);
  left.mse = immse(left.y, left.yPred);
  right.mse = immse(right.y, right.yPred);
  mse = (left.mse * numel(left.y) + right.mse * numel(right.y)) ...
    / (numel(left.y) + numel(right.y));
  gain = current.mse - mse;
end