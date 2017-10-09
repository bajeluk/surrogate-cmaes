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
  % but don't remove everything
  if isempty(trans.featuresIdx)
    trans.featuresIdx = 1:1;
  end
  X = X(:, trans.featuresIdx);
  d = size(X, 2);
  opt.nFeatures = min(opt.nFeatures, d);
  
  if opt.zscore
    [X, trans.mu, trans.sigma] = zscore(X);
  end
  featuresSampleIdx = 1:d;
  if opt.pca
    warning('off', 'stats:pca:ColRankDefX');
    [trans.coeff, X, ~, ~, explained, trans.mu] = pca(X);
    warning('on', 'stats:pca:ColRankDefX');
    % prefer features with higher explained value
    featuresSampleIdx = datasample(featuresSampleIdx, opt.nFeatures, ...
      'Replace', false, 'Weights', explained);
  else
    featuresSampleIdx = datasample(featuresSampleIdx, opt.nFeatures, ...
      'Replace', false);
  end
  trans.featuresIdx = trans.featuresIdx(featuresSampleIdx);
  valuesSampleIdx = datasample(1:n, opt.nValues, ...
    'Replace', false);
  
  X = X(valuesSampleIdx, featuresSampleIdx);
  if ~isempty(opt.polynomial) && ~strcmpi(opt.polynomial, 'linear')
    trans.polynomial = opt.polynomial;
    X = generateFeatures(X, opt.polynomial, false);
  end
  y = y(valuesSampleIdx, :);
end