function ft = feature_dispersion(X, y, settings)
% ft = FEATURE_DISPERSION(X, y, settings) returns dispersion features for 
% dataset [X, y] according to settings.
%
% The dispersion features compare the dispersion among observations within
% the initial design and among a subset of these points. The subsets are
% created based on tresholds defined using quantiles of the objective
% values. (Lunacek and Whitley, 2006)
%
% settings:
%   distance   - distance metric (similar to pdist function) | default:
%                'euclidean'
%   dist_param - additional parameter to distance (similar to pdist 
%                function)
%   minimize   - binary flag stating whether the objective function should 
%                be minimized | default: true
%   quantiles  - quantiles to calculate tresholds | default: [0.02, 0.05, 
%                0.1, 0.25]
%
% Features:
%   ratio_mean_[quantile]   - ratio of quantile and all points mean 
%                             distances
%   ratio_median_[quantile] - ratio of quantile and all points median 
%                             distances
%   diff_mean_[quantile]    - difference between subset and all points mean
%                             distances
%   diff_median_[quantile]  - difference between subset and all points 
%                             median distances

  if nargin < 3
    if nargin < 2
      help feature_dispersion
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  qnt = defopts(settings, 'quantiles', [0.02, 0.05, 0.1, 0.25]);
  min_fun = defopts(settings, 'minimize', true);
  metric = defopts(settings, 'distance', 'euclidean');
  % default dist_param value
  switch metric
    case 'mahalanobis'
      def_param = nancov(X);
    case 'minkowski'
      def_param = 2;
    case 'seuclidean'
      def_param = nanstd(X);
    otherwise
      def_param = [];
  end
  dist_param = defopts(settings, 'dist_param', def_param);
  
  nQuant = numel(qnt);
  
  if ~min_fun
    y = -y;
  end
  
  % quantile tresholds
  y_tresh = quantile(y, qnt);
  % calculate distances between points
  for q = 1:nQuant
    % all points distances
    allDist = pdist(X, metric, dist_param);
    % quantile points distances
    quantDist = pdist(X(y <= y_tresh(q), :), metric, dist_param);
    % calculate means and medians
    allDistMean = mean(allDist);
    allDistMedian = median(allDist);
    quantDistMean = mean(quantDist);
    quantDistMedian = median(quantDist);
    % calculate ratios and differences
    ft.(sprintf('ratio_mean_%02d',   100*qnt(q))) = quantDistMean   / allDistMean;
    ft.(sprintf('ratio_median_%02d', 100*qnt(q))) = quantDistMedian / allDistMedian;
    ft.(sprintf('diff_mean_%02d',    100*qnt(q))) = quantDistMean   - allDistMean;
    ft.(sprintf('diff_median_%02d',  100*qnt(q))) = quantDistMedian - allDistMedian;
  end
  
end