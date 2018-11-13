function ft = feature_cmaes(X, y, settings)
% ft = FEATURE_CMAES(X, y, settings) returns features based on CMA-ES state
% variables for dataset [X, y].
%
% settings:
%   cma_cov        - covariance matrix of CMA-ES
%   cma_evopath_c  - CMA-ES covariance matrix evolution path vector
%   cma_evopath_s  - CMA-ES step-size evolution path vector
%   cma_generation - generation number of CMA-ES algorithm
%   cma_mean       - CMA-ES mean value
%   cma_restart    - number of CMA-ES restarts
%   cma_step_size  - CMA-ES step-size value
%
% Features:

  if nargin < 3
    if nargin < 2
      help feature_cmaes
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end
  
  % parse settings
  ccov = defopts(settings, 'cma_cov', NaN);
  evopath_c = defopts(settings, 'cma_evopath_c', NaN);
  evopath_s = defopts(settings, 'cma_evopath_s', NaN);
  generation = defopts(settings, 'cma_generation', NaN);
  cmean = defopts(settings, 'cma_mean', NaN);
  restart = defopts(settings, 'cma_restart', NaN);
  step_size = defopts(settings, 'cma_step_size', NaN);
  
  dim = size(X, 2);
  
  % TODO: study CMA-ES properties and design useful features
  
  ft.cma_generation = generation;
  ft.cma_step_size = step_size;
  ft.cma_restart = restart;
  % mahalanobis distance of the CMA mean to dataset
  ft.cma_mean_dist = mahal(cmean, X);
  % determinant of covariance matrix
  ft.cma_det_cov = det(ccov);
  % norm of evolution path (used to update covariance matrix)
  %   covariance matrix p_c*p_c’ has rank one, 
  %   i.e. only one eigenvalue ||p_c||^2
  ft.cma_evopath_c_norm = norm(evopath_c);
  % step-size evolution path
  %   length of the evolution path p_\sigma compared to the expected length
  %   of random steps: ||p_\sigma|| / E (|| N(0, I) ||)
  %   E (|| N(0, I) ||) = \sqrt(2)*\gamma((D+1)/2) / \gamma(D/2)) 
  ft.cma_evopath_s_norm = norm(evopath_s) / ...
                      (sqrt(2) * gamma( (dim+1) / 2 ) / gamma( dim / 2 ));
  
end