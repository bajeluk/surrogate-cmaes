function indInv = inverseIndex(index)
  % Returns inversed index vector
  %
  % Example:
  %   inverseIndex([2, 4, 7])
  %   ans =
  %        0     1     0     2     0     0     3
  
  if nargin == 0
    help inverseIndex
    indInv = [];
    return
  end
  
  % inversed to normal
  if any(index == 0) 
    indInv = find(index);
  % normal to inversed
  else
    a = zeros(1, max(index));
    a(index) = 1;
    c = cumsum(a);
    a(index) = c(index);
    indInv = a;
  end
end