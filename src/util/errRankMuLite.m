%ERRRANKMULITE Difference in rankings between the permutations p1 and p2, counting only the 'mu' first ranked elements (of p2)
%
% err = ERRRANKMULITE(p1, p2, mu, maxErr)
%       returns the number of ordering errors between the vectors y1 and 
%       y2, but rank-errors of only the first mu ranks (ordinals according to the second
%       vectory y2) are calculated, normalized by the already calculated maxErr constant
%       p1 and p2 are expected as row-vectors
function errNorm = errRankMuLite(p1, p2, mu, maxErr)
  lambda = numel(p1);

  % speedup: the following 3 lines are the same as
  %   inRank2 = ranking(y2);
  inRank2     = zeros(1,lambda);
  inRank2(p2) = [1:lambda];

  r1 = [1:mu];
  r2 = inRank2(p1(1:mu));
  rDiff = abs(r2 - r1);

  errNorm = sum(rDiff)/maxErr;
end
