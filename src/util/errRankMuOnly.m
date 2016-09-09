function err = errRankMuOnly(x, mu)
%ERRRANKMUONLY Difference in ranks between the vector x and [1:length(x)] counting only the mu first ranked elements
%
% err = ERRRANKMUONLY(x, mu)
%       returns the number of ordering errors between the vector x and 
%       the vector [1:length(x)], but rank-errors among the first mu 
%       positions are calculated


  if (size(x,1) > 1 && size(x,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  l = length(x);
  r = ranking(x);
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

  % up to mu or floor((l-1)/2)      ---- per (l-nfirst+1) points
  nfirst = min(mu,floor((l-1)/2));
  max_err = nfirst * (l-nfirst);
  % (nfirst+1) th ... mu th         ---- per nfirst points
  max_err = max_err + max(0, mu-nfirst) * nfirst;

  % return relative ratio of errors
  err = err/max_err;
end
