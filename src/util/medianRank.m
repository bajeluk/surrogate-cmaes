function mr = medianRank(values, varargin)
% Computes ranking of vector elements, where same ranks are replaced by
% medians of shifted ranks (from preciseRank).
%
% Example:
%   a  = [ 1    5    8   13    1    8    1   21 ];
%   mr = medianRank(a);
%   mr = [ 2    4   5.5   7    2   5.5   2    8 ]

  mr = preciseRank(values, varargin{:});
  ranks = unique(mr);
  % equal ranks replace by medians
  if numel(ranks) < numel(mr)
    for r = ranks
      % average rank + actual rank - 1
      mr(mr == r) = (sum(mr == r) + 1)/2 + r - 1;
    end
  end

end