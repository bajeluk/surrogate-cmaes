function f=frosen(x)
  if size(x,1) < 2 error('dimension must be greater one'); end
  N = size(x,1); 
  popsi = size(x,2); 
  f = 1e2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + sum((x(1:end-1,:)-1).^2,1);
  % f = f + f^0.9 .* (2*randn(1,popsi) ./ randn(1,popsi).^0 / (2*N)); 
