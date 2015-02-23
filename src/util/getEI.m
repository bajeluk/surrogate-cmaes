function ei = getEI(x, y, s2, best_y)
% modelGetEI - returns Expected Improvement (EI) for a given model and input x
%
%          /    (best_y - model_y(x))*NCDF(z) + model_std(x)*NPDF(z)    if model_std(x) > 0
% EI(x) = <
%          \    0                                                       if model_std(x) == 0
% 
% ei = modelGetEI(M, x, best_y) - returns expected improvement (EI) at @x for the function 
%       modelled by the model @M, so-far reached minimum value @best_y
%
% x     vector or matrix of size n by D (n = # of samples, D = dimension)

  n = size(x,1);

  null_variance = (abs(s2) < eps);
  z = zeros(n,1);
  best_y_full = repmat(best_y,n,1);
  z(~null_variance) = (best_y_full(~null_variance) - y(~null_variance)) ./ sqrt(s2(~null_variance));

  % save the EI
  ei = zeros(n,1);
  ei(~null_variance) = (best_y_full(~null_variance) - y(~null_variance)).*normcdf(z(~null_variance)) + s2(~null_variance).*normpdf(z(~null_variance));

  function cdf = normcdf(x)
    cdf = 0.5 * erfc(-(x)/(sqrt(2)));
  end

  function pdf = normpdf(x)
    pdf = exp(-0.5 * x.^2) ./ sqrt(2*pi);
  end

end
