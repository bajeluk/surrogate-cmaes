function [X, y, trans] = transform(X, y, opt)
  [n, d] = size(X);
  opt.nFeatures = min(defopts(opt, 'nFeatures', d), d);
  opt.nValues = min(defopts(opt, 'nValues', n), n);
  opt.zscore = defopts(opt, 'zscore', false);
  opt.pca = defopts(opt, 'pca', false);
  opt.polynomial = defopts(opt, 'polynomial', 'linear');
  trans = struct;
  if opt.zscore
    [X, trans.mu, trans.sigma] = zscore(X);
  end
  if opt.pca
    [trans.coeff, X, ~, ~, explained, trans.mu] = pca(X);
    % prefer featurer with higher explained value
    trans.featuresIdx = datasample(1:d, opt.nFeatures, ...
      'Replace', false, 'Weights', explained);
  else
    trans.featuresIdx = datasample(1:d, opt.nFeatures, ...
      'Replace', false);
  end
  trans.valuesIdx = datasample(1:n, opt.nValues, ...
    'Replace', false);
  X = X(trans.valuesIdx, trans.featuresIdx);
  if ~stricmp(opt.polynomial, 'linear')
    trans.polynomial = opt.polynomial;
    X = generateFeatures(X, opt.polynomial, false);
  end
  y = y(trans.valuesIdx, :);
end