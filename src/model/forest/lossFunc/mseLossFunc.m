function err = mseLossFunc(x, y)
% ERR = MSELOSSFUNC(X,Y) calculates the mean-squared error (MSE) between 
%   the arrays X and Y.
  assert(all(size(x) == size(y)), 'mseLossFunc: Array dimensions are not equal.')
  err = (norm(x(:)-y(:), 2).^2) / numel(x);
  
end

