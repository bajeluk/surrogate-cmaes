function [XT, yt, trans] = transform(X, y, opt)
  [n, d] = size(X);
  opt.nFeatures = min(defopts(opt, 'nFeatures', d), d);
  opt.nValues = min(defopts(opt, 'nValues', n), n);
  opt.pca = defopts(opt, 'pca', false);
  opt.zscore = defopts(opt, 'zscore', false);
  trans = struct;
  XT = X;
  if opt.zscore
    [XT, trans.mu, trans.sigma] = zscore(XT);
  end
  if opt.pca
    [trans.coeff, XT, ~, ~, explained, trans.mu] = pca(XT);
    % prefer featurer with higher explained value
    trans.featuresIdx = datasample(1:d, opt.nFeatures, ...
      'Replace', false, 'Weights', explained);
  else
    trans.featuresIdx = datasample(1:d, opt.nFeatures, ...
      'Replace', false);
  end
  trans.valuesIdx = datasample(1:n, opt.nValues, ...
    'Replace', false);
  XT = XT(trans.valuesIdx, trans.featuresIdx);
  yt = y(trans.valuesIdx, :);
end