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

function X = ilhsdesign(N, dim, dup)
% ILHSDESIGN generates improved LHS design
%
% Input:
%   N   - number of points to be sampled
%   dim - number of dimensions needed
%   dup - duplication factor which affects the number of points that the
%         optimization algorithm has to choose from
%
% Beachkofski and Grandhi (2002): Improved Distributed Hypercube Sampling
% 
% Implemented according to improvedLHS_R.cpp from lhs R-package.

  if nargin < 3
    dup = 10;
  end
  
  resultMat = zeros(N, dim);
  % length of the point1 columns and the list1 vector
  len = dup*(N - 1);
  point1 = zeros(N, len);
  list1 = NaN(1, len);
  vec = NaN(1, dim);
  % optimal spacing between points
  opt = N / nthroot(N, dim);
  opt2 = opt*opt;
  
  % initialize the avail matrix
  avail = repmat(1:N, dim, 1);
  % N random numbers in the last column of result matrix
  resultMat(:, dim) = randi(N, N, 1);
  % use random numbers from the last column of resultMat to place an N
  % value randomly through the avail matrix
  avail(:, resultMat(:, N)) = N;
  
  % move backwards through the result matrix columns
  for count = N-1 : -1 : 1
    for row = 0 : dim-1
      for col = 1 : dup-1
        % create list1 vector
        for j = 0 : count - 1
          list1(j + count*col + 1) = avail(row, j);
        end
      end
      % create a set of points to choose from
      for col = count*dup : -1 : 1
        point_index = randi(col);
        point1(row + 1, col) = list1(point_index);
        list1(point_index) = list1(col);
      end
    end
    
    min_all = Inf;
    best = 0;
    for col = 0 : dup*count - 1 
      min_candidate = Inf;
      for j = count : N
        distSquared = 0;
        % find the distance between candidate points and the points already
        % in the sample
        for k = 0 : dim-1
          vec(k+1) = point1(k+1, col+1) - resultMat(k+1, j+1);
          distSquared = distSquared + vec(k+1)*vec(k+1);
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
    for row = 0 : dim-1
      resultMat(row+1, count) = point1(row+1, best);
    end
    % update the numbers that are available for the future points
    for row = 0 : dim-1
      for col = 0 : N
        if avail(row+1, col+1) == resultMat(row+1, count)
          avail(row+1, col+1) = avail(row+1, count);
        end
      end
    end
  end
  
  % once all but the last points of result are filled in, there is only one
  % choice left
  for row = 0 : dim-1
    resultMat(row+1, 1) = avail(row+1, 1);
  end

end