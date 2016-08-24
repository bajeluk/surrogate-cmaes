function err = errRank(x)
% Returns the number of errors in ranking of the vector x

  if (size(x,1) > 1 && size(x,2) > 1)
    error('Error in ranking can be done only for vectors');
  end

  [~, id] = sort(x);
  err = sum(abs(id - [1:length(x)]));
end
