function X = makeDesign(N, dim, smooth_flag, varargin)
%MAKEDESIGN Generate a design for metalearning.
%   X   = MAKEDESIGN(N, dim, smooth_flag) generate dim-dimensional sample
%   of size N. If smooth_flag is enabled, the points are placed
%   randomly within each (k*(1/N), (k+1)*(1/N)) interval.
%   X   = MAKEDESIGN(__, __, __, designType) specify design type
%   as one of 'lhs' (latin hypercube, default) or 'lhsnorm'

   if nargin > 3
     designType = varargin{1};
   else
     designType = 'lhs';
   end

   if smooth_flag
     smooth = 'on';
   else
     smooth = 'off';
   end

   if strcmp(designType, 'lhs')
     X = lhsdesign(N, dim, 'smooth', smooth);
     X = -5 + 10 * X;
   else
     mu = zeros(dim, 1);

     % scale the normal distribution such that alpha% of the mass lie
     % in the bbob's search interval [-5, 5]^dim
     alpha = 0.99;
     k = erfinv(alpha) * sqrt(2);
     sigma = (5 / k) * eye(dim);
     X = lhsnorm(mu, sigma, N, smooth);
   end
end

