function s = predictionStats(y1, y2, stat)
% s = predictionStats(y, ym, stat) count statistic for model prediction
%
% Input:
%   y1   - first data vector
%   y2   - second data vector
%   stat - statistic to compute:
%            kendall  - Kendall' tau rank
%            mse      - mean square error
%            mzoe     - mean zero-one error
%            rankmse  - mean square error of ranks
%            rankmzoe - mean zero-one error of ranks
%            rde      - ranking difference error
%            mae      - mean absolute error
%            r2       - r squared, aka explained variance, aka coeff. of
%                       determination

  % normalize input
  y1 = y1(:);
  y2 = y2(:);

  n = numel(y1);
  assert(numel(y2) == n, 'Lengths of compared vectors has to be equal.')

  % compute ranking statistics
  if length(stat) > 3 && strcmp(stat(1:4), 'rank')
    [~, id] = sort(y1);
    y1(id) = (1:n);
    [~, id] = sort(y2);
    y2(id) = (1:n);
    stat = stat(5:end);
  end

  switch stat
    case 'mse'
      s = sum((y2 - y1).^2) / n;
    case 'mzoe'
      s = sum(y2 ~= y1) / n;
    case 'kendall'
      s = corr(y2, y1, 'type', 'kendall');
    case 'rde'
      s = errRankMu(y2, y1, floor(n/2));
    case 'mae'
      s = sum(abs(y2 - y1)) / n;
    case 'r2'
      y1_mean = mean(y1);
      % regression sum of squares
      ss_reg = sum((y2 - y1_mean).^2);
      % total sum of squares
      ss_tot = sum((y1 - y1_mean).^2);
      s = ss_reg / ss_tot;
    otherwise
      s = NaN;
  end

end