function tr = tolerantRank(values, tol)
% Computes ranking of vector elements, where number of equal ranks does not
% play role.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   tr = tolerantRank(a);
%   tr = [ 1    2    3    4    1    3    1    5 ]

  if (exist('tol', 'var') && ~isempty(tol))
    [~, ~, tr] = uniquetol(values, tol);
  else
    [~, ~, tr] = unique(values);
  end
  tr = tr';

end