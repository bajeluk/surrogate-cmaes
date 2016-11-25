%ERRRANKMULITE Difference in rankings between permutations si1 and si2, counting only the 'mu' first ranked elements of si2
%
% err = ERRRANKMULITE(si1, si2, mu, maxErr)
%       returns the number of ordering errors between the permutations si1 and 
%       si2, but rank-errors of only the first mu ranks (ordinals according to the second
%       permutation si2) are calculated, normalized by the already calculated maxErr constant
%       si1 and si2 are expected as row-vector results of [~, si] = sort(...)
function errNorm = errRankMuLite(si1, si2, mu, maxErr)
  lambda = numel(si1);

  % speedup: the following 2 lines are the same as
  %   inRank2 = ranking(si2);
  inRank1     = zeros(1,lambda);
  inRank1(si1) = (1:lambda);

  r1 = inRank1(si2(1:mu));
  r2 = (1:mu);
  errNorm = sum(abs(r2 - r1))/maxErr;
end
