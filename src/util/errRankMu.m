function err = errRankMu(x, mu)
%ERRRANKMU Difference in ranks between the vector x and [1:length(x)] counting the (lambda-mu) last elements only for 1 errorpoint
%
% err = ERRRANKMU(x, mu)
%       returns the number of ordering errors between the vector x and 
%       the vector [1:length(x)], but the ordering among elements
%       (mu+1),(mu+2)...(length(x)) is not considered as being error
%       and positioning them among the first mu positions is calculated
%       for only one error-point


  if (size(x,1) > 1 && size(x,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  l = length(x);
  r = ranking(x);
  r = r(:);
  % calculate the differences to the right ranking
  err = abs(r - [1:l]');
  % those errors which belong to ranking > mu calculate only
  % for 1 point (instead of the difference to the right ranking)
  err((err > 0) & (r > mu)) = 1;
  % do not calculate errors in the order between elements of rank
  % higher than mu that are actually in position higher than mu
  mu_mask = false(l,1);
  mu_mask(mu+1:end) = true;
  err((r > mu) & mu_mask) = 0;
  % return the overall number of errors
  err = sum(err);

  % Try to calculate maximum number of possible errors
  %
  % but it is not totally right, at least for the case l=5, mu=4
  %
  % up to mu or floor((l-1)/2)      ---- per (l-nfirst+1) points
  nfirst = min(mu,floor((l-1)/2));
  max_err = nfirst * (l-nfirst+1);
  % (nfirst+1) th ... mu th         ---- per nfirst points
  max_err = max_err + max(0, mu-nfirst) * nfirst;

  % return relative ratio of errors
  err = err/max_err;
end
