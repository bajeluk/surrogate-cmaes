function d = kendall_dist_mu(x, y, mu)
%KENDALL Kendall's tau distance between two vectors x and y
%
% err = KENDALL_DIST(x, y)
%       returns Kendall's tau distance between the vector x 
%       and the vector y


  if (size(x,1) > 1 && size(x,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  % calculates Kendall's tau distance:
  l = length(x);
  d = 0;
  for i = 1:(l-1)
    for j = (i+1):l
      d = d + ((x(i) < x(j)) == ~(y(i) < y(j)));
    end
  end
end
