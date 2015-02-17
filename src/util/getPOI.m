function poi = getPOI(x, y, s2, target)
% getPOI - returns the probability of improvement (POI)
% 
% poi = getPOI(y, s2, x, target)
%   returns probability of improvement (POI) for specified @target (usually 
%   somewhat better than the best f(x) sofar)
%   at @x's where the model predicted @y with variance @s2
% 
% x     vector or matrix of size n by D (n = # of samples, D = dimension)

  dim = size(x, 2);
  n = size(x,1);

  % evaluate the PoI
  % be carefull when dividing by 0 (edit 17/10/2013)
  null_variance = (abs(s2) < eps);
  q = zeros(n,1);
  target_full = repmat(target,n,1);
  q(~null_variance) = (target_full(~null_variance) - y(~null_variance)) ./ sqrt(s2(~null_variance));

  % save the PoI
  poi = zeros(n,1);
  poi(~null_variance) = normcdf(q(~null_variance));

  function cdf = normcdf(x)
    cdf = 0.5 * erfc(-(x)/(sqrt(2)));
  end

end
