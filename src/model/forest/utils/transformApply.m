function [X] = transformApply(X, trans)
  if isfield(trans, 'mu')
    X = bsxfun(@minus, X, trans.mu);
  end
  if isfield(trans, 'sigma')
    X = bsxfun(@rdivide, X, trans.sigma);
  end
  if isfield(trans, 'coeff')
    X = X * trans.coeff;
  end
  X = X(:, trans.featuresIdx);
  if isfield(trans, 'polynomial')
    X = generateFeatures(X, trans.polynomial, false);
  end
end