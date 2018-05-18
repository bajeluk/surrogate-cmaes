function s = logsumexp(x, varargin)
%LOGSUMEXP Computes log(sum(w .* exp(x), dim)). Useful for computing sum of
%densities.
%   x   -- an array or a matrix
%   dim -- dimension to operate on
%   w   -- vector of weights

  if nargin > 2
    dim = varargin{2};
  else
    if size(x, 1) == 1
      dim = 2;
    else
      dim = 1;
    end
  end

  if nargin > 1
    w = varargin{1};
    if isempty(w)
      w = ones(size(x));
    elseif (ndims(w) > 1 && ~all(size(w) == size(x))) || ...
        (ndims(w) == 1 && ~any(size(w) == size(x)))
      error('logsumexp: Dimensions of the input and the weights must agree.');
    end
  else
    w = ones(size(x));
  end

  [m, i] = max(x, [], dim);
  m(isinf(m)) = 0;

  % relies of automatic dimension expansion of '-' and '.*' operators
  s = log(sum(w .* exp(x - m), dim)) + m;

end

