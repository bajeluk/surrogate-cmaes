function mr = medianRank(values)
% Computes ranking of vector elements, where same ranks are replaced by
% medians of shifted ranks (from preciseRank).
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   mr = medianRank(a);
%   mr = [ 2    4   5.5   7    2   5.5   2    8 ]

  mr = preciseRank(values);
  ranks = unique(mr);
  % equal ranks replace by medians
  if numel(ranks) < numel(mr)
    for r = ranks
      % average rank + actual rank - 1
      mr(mr == r) = (sum(mr == r) + 1)/2 + r - 1;
    end
  end

end

function tr = tolerantRank(values)
% Computes ranking of vector elements, where number of equal ranks does not
% play role.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   tr = tolerantRank(a);
%   tr = [ 1    2    3    4    1    3    1    5 ]

  [~, ~, tr] = unique(values);
  tr = tr';

end

function pr = preciseRank(values)
% Computes ranking of vector elements, where number of equal ranks shifts
% the following ranks.
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   pr = preciseRank(a);
%   pr = [ 1    4    5    7    1    5    1    8 ]

  tr = tolerantRank(values);
  pr = arrayfun(@(x) sum(tr < x), tr) + 1;

end