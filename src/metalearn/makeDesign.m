function X = makeDesign(N, dim, varargin)
%MAKEDESIGN Generate a design for metalearning.
%   X   = MAKEDESIGN(dim, N) generate dim-dimensional sample
%   of size N.
%   X   = MAKEDESIGN(__, __, designType) specify design type
%   as one of 'lhs' (latin hypercube, default) or 'lhsnorm'

   if nargin > 2
     designType = varargin{1};
   else
     designType = 'lhs';
   end

   if strcmp(designType, 'lhs')
     X = lhsdesign(N, dim);
   else
     mu = zeros(dim, 1);
     sigma = eye(dim);
     X = lhsnorm(mu, sigma, N, 'on');
   end
end

