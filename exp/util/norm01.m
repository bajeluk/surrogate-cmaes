%NORM01 Normalize vector to 0,1
function y = norm01(x)
  y = (x - min(x)) ./ (max(x) - min(x));
end
