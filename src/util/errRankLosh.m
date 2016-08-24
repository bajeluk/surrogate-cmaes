function err = errRankLosh(x)
% Returns the number of errors in ranking of the vector x

  if (size(x,1) > 1 && size(x,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  l = length(x);
  err = 0;
  for i = 1:l
    for j = (i+1):l
      err = err + (x(i) > x(j));
    end
  end
  err = - (err / (l*(l-1)) * 4 - 1);
end
