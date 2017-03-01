function o = ranking(x)
% Returns ranking of the elements of the vector x  

  if (size(x,1) > 1 && size(x,2) > 1)
    error('Ranking can be done only for vectors');
  end

  [~, id] = sort(x);
  o = zeros(size(x));
  o(id) = 1:length(x);
end
