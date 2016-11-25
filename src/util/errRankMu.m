%ERRRANKMU Difference in rankings between the vectors y1 and y2, counting only the 'mu' first ranked elements (of y2)
%
% err = ERRRANKMU(y1, y2, mu)
%       returns the number of ordering errors between the vectors y1 and 
%       y2, but rank-errors of only the first mu ranks (ordinals according to the second
%       vectory y2) are calculated
%
function [errNorm, errSum, maxErr] = errRankMu(y1, y2, mu)

  % for saving already computed errRankMu normalizing constants
  persistent maxErrs;

  % if ((size(y1,1) > 1 && size(y1,2) > 1) || (size(y2,1) > 1 && size(y2,2) > 1) ...
  %     || any(size(y1) ~= size(y2)))
  %   error('Error in ranking can be done only for two vectors of same size');
  % end

  y1 = y1(:);
  y2 = y2(:);
  if (numel(y1) ~= numel(y2))
    error('Error in ranking can be done only for two vectors of same size');
  end

  if (mu == 0)
    errSum = 0;
    warning('Calling errRankMu() with mu==0 does not really make sense.');
    return;
  end

  lambda = length(y1);
  [~, si_y1] = sort(y1);

  % speedup: the following 3 lines are the same as
  %   inRank2 = ranking(y2);
  [~, si_y2] = sort(y2);
  inRank2 = zeros(lambda,1);
  inRank2(si_y2) = [1:lambda]';

  r1 = [1:lambda]';
  r2 = inRank2(si_y1);
  rDiff = abs(r2 - r1);
  rDiff((mu+1):end) = 0;

  errSum = sum(rDiff);

  %
  % Normalize the error to the range [0,1]
  %

  if (~isempty(maxErrs) && length(maxErrs) >= lambda && maxErrs(lambda) > 0)
    % we have the normalizing constant already computed
    maxErr = maxErrs(lambda);
  else
    if (isempty(maxErrs))
      maxErrs = zeros(1,lambda);
    end
    % First, calculate maximal Ranking differences to the ranking 1:n
    % for every point for two cases, i.e. when
    % 1) going into opposite end of the not-so-far occupied part
    going_opposite = abs(-2*[1:lambda] + lambda+1);
    going_opposite((mu+1):end) = 0;
    % 2) going to the left (also not further than not-yet occupied part)
    maxErr = 0;
    for ri = 2:lambda
      % we have to iterate through possible positions of split between
      % 'opposite' and 'left' part
      going_left = min([1:lambda] - 1, ri - 1);
      % Maximal error for each point in the first part would be going opposite
      max_err = going_opposite;
      % But the second part, starting from the reverting index,
      % is not clear: decide which option yields more error points
      if (sum(going_opposite(ri:end)) < sum(going_left(ri:end)))
        max_err(ri:end) = going_left(ri:end);
      end
      % again, calculate only the cases where y1 <= mu
      max_err((mu+1):end) = 0;
      % Sum all the errors
      sum_max_err = sum(max_err);
      if (sum_max_err > maxErr)
        maxErr = sum_max_err;
      end
    end
    maxErrs(lambda) = maxErr;
  end

  % return relative ratio of errors
  errNorm = errSum/maxErr;
end
