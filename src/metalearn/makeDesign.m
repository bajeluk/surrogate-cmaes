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

   if nargin < 3 || smooth_flag
     smooth = 'on';
   else
     smooth = 'off';
   end
   
   % half of bbob's search interval [-5, 5]^dim in one dimension
   int_range = 5;

   if strcmp(designType, 'lhs')
     X = lhsdesign(N, dim, 'smooth', smooth);
     X = -int_range + 2 * int_range * X;
   elseif strcmp(designType, 'ilhs')
     X = ilhsdesign(N, dim, 1, smooth);
     X = -int_range + 2 * int_range * X;
   else
     mu = zeros(dim, 1);

     % scale the normal distribution such that alpha% of the mass lie
     % in the bbob's search interval [-5, 5]^dim
     alpha = 0.99;
     k = erfinv(alpha) * sqrt(2);
     sigma = (int_range / k) * eye(dim);
     X = lhsnorm(mu, sigma, N, smooth);
   end
end

function X = ilhsdesign(N, dim, dup, smooth)
% ILHSDESIGN generates improved LHS design
%
% Input:
%   N      - number of points to be sampled
%   dim    - number of dimensions needed
%   dup    - duplication factor which affects the number of points that the
%            optimization algorithm has to choose from | {1,2,...}
%   smooth - 'off' produces points at the midpoints of the intervals:
%            0.5/n, 1.5/n, ..., 1-0.5/n | {'on', 'off'}
%
% Beachkofski and Grandhi (2002): Improved Distributed Hypercube Sampling
% 
% Implemented according to improvedLHS_R.cpp from lhs R-package.

  if nargin < 4
    if nargin < 3
      dup = 1;
    end
    smooth = 'on';
  end
  assert(dup == floor(dup) && dup > 0 && isfinite(dup), ...
         'dup has to be a finite positive integer')
  
  lhs_mat = zeros(dim, N);
  % length of the point1 columns and the list1 vector
  len = dup*(N - 1);
  point1 = NaN(dim, len);
  list1 = NaN(1, len);
  vec = NaN(1, dim);
  % optimal spacing between points
  opt = N / nthroot(N, dim);
  opt2 = opt*opt;
  
  % initialize the avail matrix
  avail = repmat(1:N, dim, 1);
  % N random numbers in the last column of result matrix
  lhs_mat(:, N) = randi(N, dim, 1);
  % use random numbers from the last column of lhs_mat to place an N
  % value randomly through the avail matrix
  for row = 1:dim
    avail(row, lhs_mat(row, N)) = N;
  end
  
  % move backwards through the result matrix columns
  for count = N-1 : -1 : 1
    for row = 1 : dim
      for col = 0 : dup-1
        % create list1 vector
        for j = 1 : count
          list1(j + count*col) = avail(row, j);
        end
      end
      % create a set of points to choose from
      for col = count*dup : -1 : 1
        point_index = randi(col);
        point1(row, col) = list1(point_index);
        list1(point_index) = list1(col);
      end
    end
    
    min_all = Inf;
    best = 1;
    for col = 1 : dup*count - 1
      min_candidate = Inf;
      for j = count+1 : N
        distSquared = 0;
        % find the distance between candidate points and the points already
        % in the sample
        for k = 1 : dim
          vec(k) = point1(k, col) - lhs_mat(k, j);
          distSquared = distSquared + vec(k)*vec(k);
        end
        % if distSquared is the smallest so far replace it in min_candidate
        if min_candidate > distSquared
          min_candidate = distSquared;
        end
      end
      % if the difference between min_candidate and opt2 is the smallest so
      % far, then keep that point as the best
      if abs(min_candidate - opt2) < min_all
        min_all = abs(min_candidate - opt2);
        best = col;
      end
    end
    
    % take the best point out of point1 and place it in the result matrix
    lhs_mat(:, count) = point1(:, best);
    % update the numbers that are available for the future points
    for row = 1 : dim
      for col = 1 : N
        if avail(row, col) == lhs_mat(row, count)
          avail(row, col) = avail(row, count);
        end
      end
    end
  end
  
  % once all but the last points of result are filled in, there is only one
  % choice left
  lhs_mat(:, 1) = avail(:, 1);
  
  % generate exact points from lhs_mat according to smoothing strategy
  lhs_int = 1/N;
  if strcmp(smooth, 'on')
    X = lhs_int*(lhs_mat' - 1 + rand(N, dim));
  else
    X = lhs_int*(lhs_mat' - 0.5);
  end

end