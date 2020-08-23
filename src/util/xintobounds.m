function [x, idx] = xintobounds(x, lbounds, ubounds)
% XINTOBOUNDS replaces elements of column vector x (or matrix consisting of 
% column vectors) out of lower and upper bounds (lbounds, ubounds) by
% values from lbounds, ubounds.
%
% [x, idx] = XINTOBOUNDS(x, lbounds, ubounds)
%
% Input:
%   x       - column vector or a matrix consisting of column vectors |
%             double
%   lbounds - vector of lower bounds | double
%   ubounds - vector of upper bounds | double
%
% Output:
%   x   - column or a matrix of columns respecting bounds | double
%   idx - rebounding indices | integer vector
%
% See Also:
%   cmaes, s_cmaes

  if nargin < 3
    help xintobounds
    return
  end

  if ~isempty(lbounds)
    if length(lbounds) == 1
      idx = x < lbounds;
      x(idx) = lbounds;
    else
      arbounds = repmat(lbounds, 1, size(x,2));
      idx = x < arbounds;
      x(idx) = arbounds(idx);
    end
  else
    idx = 0;
  end
  if ~isempty(ubounds)
    if length(ubounds) == 1
      idx2 = x > ubounds;
      x(idx2) = ubounds;
    else
      arbounds = repmat(ubounds, 1, size(x,2));
      idx2 = x > arbounds;
      x(idx2) = arbounds(idx2);
    end
  else
    idx2 = 0;
  end
  idx = idx2-idx; 
  
end