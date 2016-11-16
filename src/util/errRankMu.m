%ERRRANKMU Difference in rankings between the vectors y1 and y2, counting only the 'mu' first ranked elements (of y2)
%
% err = ERRRANKMU(y1, y2, mu)
%       returns the number of ordering errors between the vectors y1 and 
%       y2, but rank-errors of only the first mu ranks (ordinals according to the second
%       vectory y2) are calculated
%
function err = errRankMu(y1, y2, mu)
  if ((size(y1,1) > 1 && size(y1,2) > 1) || (size(y2,1) > 1 && size(y2,2) > 1) ...
      || any(size(y1) ~= size(y2)))
    error('Error in ranking can be done only for two vectors of same size');
  end

  if (mu == 0)
    err = 0;
    warning('Calling errRankMu() with mu==0 does not really make sense.');
    return;
  end

  lambda = length(y1);
  [~, si1] = sort(y1);
  [~, si2] = sort(y2);
  rDiff = abs(si2 - si1);
  rDiff((mu+1):end) = 0;

  err = sum(rDiff);

  %
  % Normalize the error to the range [0,1]
  %
  
  % First, calculate maximal Ranking differences to the ranking 1:n
  % for every point for two cases, i.e. when
  % 1) going into oposite end of the not-so-far occupied part
  going_oposite = abs(-2*[1:lambda] + lambda+1);
  going_oposite((mu+1):end) = 0;
  % 2) going to the left (also not further than not-yet occupied part)
  riMaxErr = 0;
  for ri = 2:lambda
    going_left = min([1:lambda] - 1, ri - 1);
    % Maximal error for each point in the first part would be going oposite
    max_err = going_oposite;
    % But the second part, starting from the reverting index,
    % is not clear: decide which option yields more error points
    if (sum(going_oposite(ri:end)) < sum(going_left(ri:end)))
      max_err(ri:end) = going_left(ri:end);
    end
    % again, calculate only the cases where y1 <= mu
    max_err((mu+1):end) = 0;
    % Sum all the errors
    sum_max_err = sum(max_err);
    if (sum_max_err > riMaxErr)
      riMaxErr = sum_max_err;
    end
  end

  % return relative ratio of errors
  err = err/riMaxErr;
end
