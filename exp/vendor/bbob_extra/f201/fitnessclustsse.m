function fval = fitnessclustsse(centers, data)
% fval = fitnessclustsse(centers, data)
%
% Function fitnessclustsse implements the Euclidean sum of squares 
% objective function for the clustering problems.
%
% Input: 
%   centers - (kp x 1) candidate solution vector representing the positions
%             of k cluster centers. That is, the dimensionality of the 
%             optimization problem is kp.
%   data    - (n x p) dataset to be clustered (n data points of 
%             dimensionality p)
%
% Output:
%   fval - objective function value
%
% Note: This is the upgraded version of the original from 
% http://realopt.uqcloud.net/ess_clustering.html

  % parse input
  if nargout > 0
    fval = NaN;
  end
  if nargin < 2
    help fitnessclustsse
    return
  end
  
  [kp, m] = size(centers);
  [n, p]  = size(data);
  assert(mod(kp, p) == 0, ...
    'Candidate solution vector lenght differs from k*%d, where k is natural', p);
  k = kp/p;

  % Need to convert "centers" from kp x 1 candidate solution vector to a 
  % center coordinates - k x p
  centers = reshape(centers, p, k);
  centers = centers';

  % Pairwise distances of all data points
  D = pdist2(centers, data);

  % Assign closest cluster centers to each data point
  [~, ind] = min(D, [], 1);

  % Calculate error value for current clustering
  errorclust = 0;
  for i=1:k
      errorclust = errorclust + sum( pdist2(centers(i,:),data(ind==i,:)).^2 );
  end

  fval = errorclust;
end