function [XT] = transformApply(X, trans)
  XT = X;
  if isfield(trans, 'mu')
    XT = bsxfun(@minus, XT, trans.mu);
  end
  if isfield(trans, 'sigma')
    XT = bsxfun(@rdivide, XT, trans.sigma);
  end
  if isfield(trans, 'coeff')
    XT = XT * trans.coeff;
  end
  XT = XT(:, trans.featuresIdx);
end