function pr = preciseRank(values, varargin)
% Computes ranking of vector elements, where number of equal ranks shifts
% the following ranks.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   pr = preciseRank(a);
%   pr = [ 1    4    5    7    1    5    1    8 ]

  tr = tolerantRank(values, varargin{:});
  pr = arrayfun(@(x) sum(tr < x), tr) + 1;

end