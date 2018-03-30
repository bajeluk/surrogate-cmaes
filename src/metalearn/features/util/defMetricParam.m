function def_param = defMetricParam(metric, X)
% defMetricParam(metric) returns default values for distance metric 
% parameters according to pdist2.
%
% Input:
%   metric - pdist2 distance parameter
%   X      - pdist2 input data (necessary for some metrics)
%
% Output:
%
%       metric        parameter
%   ----------------------------
%    'mahalanobis'    nancov(X)
%    'minkowski'         2
%    'seuclidean'     nanstd(X)
%
% See Also:
%   pdist2

  if nargin < 2
    X = [];
  end
  
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
end