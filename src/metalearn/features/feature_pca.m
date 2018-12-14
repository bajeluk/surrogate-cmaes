function ft = feature_pca(X, y, settings)
% ft = FEATURE_PCA(X, y, settings) returns principal component features for 
% dataset [X, y] according to settings.
%
% The features extract information from a principal component analysis, 
% which is performed on the initial design – both, including and excluding 
% the objective values – and that are based on the covariance, as well as 
% correlation matrix, respectively. For each of these four combinations, 
% the features measure the proportion of principal components that are 
% needed to explain a pre-defined percentage of the data’s variance, as 
% well as the percentage of variation that is explained by the first 
% principal component. (Kershke, 2017)

%
% settings:
%   cov_x     - percentage of the data's variance to explain covariance
%               matrix of X
%   corr_x    - percentage of the data's variance to explain correlation
%               matrix of X
%   cov_init  - percentage of the data's variance to explain covariance
%               matrix of [X, y]
%   corr_init - percentage of the data's variance to explain correlation
%               matrix of [X, y]
%
% Features:
%   pca_cov_x         - proportion of principal components that are needed 
%                       to explain a pre-defined percentage of the X’s 
%                       variance
%   pca_corr_x        - proportion of principal components that are needed 
%                       to explain a pre-defined percentage of the X’s 
%                       correlation
%   pca_cov_init      - proportion of principal components that are needed 
%                       to explain a pre-defined percentage of the [X, y]’s 
%                       variance
%   pca_corr_init     - proportion of principal components that are needed 
%                       to explain a pre-defined percentage of the [X, y]’s 
%                       correlation
%   pca_pc1_cov_x     - percentage of variation that is explained by the 
%                       first principal component using covariance of X
%   pca_pc1_corr_x    - percentage of variation that is explained by the 
%                       first principal component using correlation of X
%   pca_pc1_cov_init  - percentage of variation that is explained by the 
%                       first principal component using covariance of 
%                       [X, y]
%   pca_pc1_corr_init - percentage of variation that is explained by the 
%                       first principal component using correlation of 
%                       [X, y]

  if nargin < 3
    if nargin < 2
      help feature_pca
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  settings.cov_x = defopts(settings, 'cov_x', 0.9);
  settings.corr_x = defopts(settings, 'corr_x', 0.9);
  settings.cov_init = defopts(settings, 'cov_init', 0.9);
  settings.corr_init = defopts(settings, 'corr_init', 0.9);
  
  [nData, dim] = size(X);
  
  % initial design with the objective values
  X_init = [X, y];
  
  % calculate explaining variances
  if isempty(X)
    cov_x = NaN;
    corr_x = NaN;
  else
    cov_x = explainVariance(X, @cov);
    if nData < 2
      corr_x = NaN;
    else
      corr_x = explainVariance(X, @corr);
    end
  end
  
  % calculate features
  ft.pca_cov_x = find(cov_x >= settings.cov_x, 1, 'first') / dim;
  ft.pca_corr_x = find(corr_x >= settings.corr_x, 1, 'first') / dim;
  ft.pca_pc1_cov_x = cov_x(1);
  ft.pca_pc1_corr_x = corr_x(1);
  % calculate y-dependent features
  if isempty(y) || any(isnan(y))
    ft.pca_cov_init = NaN;
    ft.pca_corr_init = NaN;
    ft.pca_pc1_cov_init = NaN;
    ft.pca_pc1_corr_init = NaN;
  else
    % calculate explaining variances
    cov_init = explainVariance(X_init, @cov);
    if nData > 1
      corr_init = explainVariance(X_init, @corr);
    else
      corr_init = NaN;
    end
    % calculate features
    ft.pca_cov_init = find(cov_init >= settings.cov_init, 1, 'first') / (dim + 1);
    ft.pca_corr_init = find(corr_init >= settings.corr_init, 1, 'first') / (dim + 1);
    ft.pca_pc1_cov_init = cov_init(1);
    ft.pca_pc1_corr_init = corr_init(1);
  end
  
  % ensure features to be non-empty in case of low number of points
  if nData < 2
    ft = repStructVal(ft, @isempty, NaN, 'test');
  end
  
end

function ev = explainVariance(X, fun)
% calculate explaining fun
  ev = fliplr(eig(fun(X)));
  ev = cumsum(ev) / sum(ev);
end