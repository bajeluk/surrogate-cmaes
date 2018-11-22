function ft = feature_nearest_better(X, y, settings)
% ft = FEATURE_NEAREST_BETTER(X, y, settings) returns nearest better
% clustering features for dataset [X, y] according to settings.
%
% Nearest better clustering features extract information based on the 
% comparison of the sets of distances from all observations towards their 
% nearest neighbors (D_nn) and their nearest better neighbors (D_nb) (see 
% Kerschke et al. 2015).
%
% The distance to the nearest neighbor of a search point x from a 
% population P:
%     d_nn(x, P) = min({dist(x, y) | y \in P \ {x}})
%
% The distance to the nearest better neighbor (f(x) stands for the 
% objective function value that corresponds to x):
%     d_nb(x, P) = min({dist(x, y) | (f(y) < f(x)) AND (y \in P)})
%
% The set of all nearest neighbor distances within the population:
%     D_nn = {d_nn(x, P) | x \in P}
%
% The set of all nearest-better distances:
%     D_nb = {d_nb(x, P) | x \in P}
%
% settings:
%   distance   - distance metric (similar to pdist function) | default:
%                'euclidean'
%   dist_param - additional parameter to distance (similar to pdist 
%                function)
%   minimize   - binary flag stating whether the objective function should 
%                be minimized | default: true
%
% Features:
%   nb_std_ratio   - ratio of the standard deviations between the D_nn and
%                    D_nb
%   nb_mean_ratio  - ratio of the means between the D_nn and D_nb
%   nb_cor         - correlation between the distances of the nearest 
%                    neighbors and nearest better neighbors
%   dist_ratio     - coefficient of variation of the distance ratios
%   nb_fitness_cor - correlation between the fitness value, and the count 
%                    of observations to whom the current observation is the
%                    nearest better neighbour, i.e., the so-called 
%                    “indegree”

  if nargin < 3
    if nargin < 2
      help feature_nearest_better
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
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
  
  % calculate distances
  distMat = squareform(pdist(X, metric, dist_param));
  % replace diagonal zeros by NaN values
  distMat = distMat + diag(NaN(1, numel(y)));
  % set of all distances to the nearest neighbors
  dnn = min(distMat)';  
  % nearest better neighbor
  [dnb, nb_indegree] = nearestBetterNeighbor(distMat, (2*min_fun - 1) * y, dnn);
  % calculate quotient of an individual's nearest-neighbor and
  % nearest-better distances
  qnn_nb = dnn/dnb;

  % calculate features
  ft.nb_std_ratio = std(dnn) / std(dnb);
  ft.nb_mean_ratio = mean(dnn) / mean(dnb);
  ft.nb_cor = corr(dnn, dnb);
  ft.dist_ratio = std(qnn_nb) / mean(qnn_nb);
  ft.nb_fitness_cor = corr(y, nb_indegree);
  
  % ensure features to be non-empty in case of empty input
  if isempty(X) || isempty(y)
    ft = repStructVal(ft, @isempty, NaN, 'test');
  end
  
end

function [dnb, nb_indegree] = nearestBetterNeighbor(distMat, y, dnn)
% calculate nearest better neighbor distances and indegrees

  nData = numel(y);
  dnb = NaN(nData, 1);
  nb_indegree = zeros(nData, 1);
  
  for i = 1:nData
    % find minimal distances
    [m, ~, id] = unique(distMat(i, y(i) > y));
    if ~isempty(m)
      % set of all distances to the nearest better neighbor
      dnb(i) = m(1);    
      % number of search points for which a certain search point is the 
      % nearest better point
      nb_indegree(i) = sum(id == 1);
    end
  end
  % replace distances for points with minimal y by dnn values
  dnb(y == min(y)) = dnn(y == min(y));
end