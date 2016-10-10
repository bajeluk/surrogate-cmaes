%ERRRANKMU Difference in rankings between the vectors y1 and y2, counting only the 'mu' first ranked elements (of y2)
%
% err = ERRRANKMU(y1, y2, mu)
%       returns the number of ordering errors between the vectors y1 and 
%       y2, but rank-errors of only the first mu ranks (ordinals according to the second
%       vectory y2) are calculated
%
function err = errRankMu(y1, y2, mu)
  [~, sort1] = sort(y1);
  ranking2   = ranking(y2);
  r = ranking2(sort1);

  if (size(r,1) > 1 && size(r,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  l = length(r);
  r = r(:);
  % calculate the differences to the right ranking
  err = abs(r - [1:l]');
  % do not calculate errors which belong to ranking > mu
  err((err > 0) & (r > mu)) = 0;
  % return the overall number of errors
  err = sum(err);

  % Try to calculate maximum error value
  %
  % it should be situation (for nfirst = min(mu,floor((l-1)/2))  )
  %
  % [(nfirst+1) (nfirst+2) ... (nfirst) ... 3 2 1]
  %
  % but it it not 100% true, it happend that 1.01 error was returned

  % up to mu or floor((l-1)/2)      ---- per (l-nfirst+1) points
  nfirst = min(mu,floor((l-1)/2));
  max_err = nfirst * (l-nfirst);
  % (nfirst+1) th ... mu th         ---- per nfirst points
  max_err = max_err + max(0, mu-nfirst) * nfirst;

  % return (approximate) relative ratio of errors
  err = err/max_err;
end
