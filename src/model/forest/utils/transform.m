function [X, y, trans] = transform(X, y, opt)
  [n, d] = size(X);
  opt.nFeatures = min(defopts(opt, 'nFeatures', d), d);
  opt.nValues = min(defopts(opt, 'nValues', n), n);
  opt.zscore = defopts(opt, 'zscore', false);
  opt.pca = defopts(opt, 'pca', false);
  opt.polynomial = defopts(opt, 'polynomial', '');
  trans = struct;
  
  % remove constant features
  varX = var(X);
  trans.featuresIdx = find(varX >= eps(max(varX)) * n);
  X = X(:, trans.featuresIdx);
  d = size(X, 2);
  opt.nFeatures = min(opt.nFeatures, d);
  
  if opt.zscore
    [X, trans.mu, trans.sigma] = zscore(X);
  end
  if opt.pca
    warning('off', 'stats:pca:ColRankDefX');
    [trans.coeff, X, ~, ~, explained, trans.mu] = pca(X);
    warning('on', 'stats:pca:ColRankDefX');
    % prefer features with higher explained value
    trans.featuresIdx = datasample(trans.featuresIdx, opt.nFeatures, ...
      'Replace', false, 'Weights', explained);
  else
    trans.featuresIdx = datasample(trans.featuresIdx, opt.nFeatures, ...
      'Replace', false);
  end
  valuesIdx = datasample(1:n, opt.nValues, ...
    'Replace', false);
  X = X(valuesIdx, trans.featuresIdx);
  if ~isempty(opt.polynomial) && ~strcmpi(opt.polynomial, 'linear')
    trans.polynomial = opt.polynomial;
    X = generateFeatures(X, opt.polynomial, false);
  end
  y = y(valuesIdx, :);
end