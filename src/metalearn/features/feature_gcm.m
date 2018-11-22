function ft = feature_gcm(X, y, settings)
% ft = FEATURE_GCM(X, y, settings) returns generalized cell mapping 
% features for dataset [X, y] according to settings.
%
% The generalized cell mapping features discretizes the continuous input 
% space utilizing a pre-defined number of blocks (cells) per dimension. 
% Each cell is represented by exactly one of its observations (best, 
% average value, or the one located closest to the cell center) and each of
% the cells is then considered to be an absorbing Markov chain. That is, 
% for each cell the transition probability for moving from one cell to one 
% of its neighboring cells is computed. Based on the resulting transition
% probabilities, the cells are categorized into attractor, transient, 
% periodic and uncertain cells.
% (Kerschke et al., 2014)
%
% settings:
%   approach   - type of representative observation | {'min', 'mean',
%                'near'}
%   blocks     - number of cell blocks per dimension
%   distance   - distance metric (similar to pdist function) | default:
%                'euclidean'
%   dist_param - additional parameter to distance (similar to pdist 
%                function)
%   lb         - lower bounds of the input space
%   minimize   - binary flag stating whether the objective function should 
%   ub         - upper bounds of the input space
%
% Features:
%   [approach]_attractors             - total number of attractor cells 
%   [approach]_pcells                 - ratio of periodic cells
%   [approach]_tcells                 - ratio of transient cells
%   [approach]_uncertain              - ratio of uncertain cells
%   [approach]_basin_prob_min         - minimum of the probabilities for
%                                       reaching the different basins
%   [approach]_basin_prob_mean        - mean of the probabilities for
%                                       reaching the different basins
%   [approach]_basin_prob_median      - median of the probabilities for
%                                       reaching the different basins
%   [approach]_basin_prob_max         - maximum of the probabilities for
%                                       reaching the different basins
%   [approach]_basin_prob_std         - standard deviation of the 
%                                     - probabilities for reaching the 
%                                       different basins
%   [approach]_basin_certain_min      - minimum of the basin sizes, basins
%                                       formed by the "certain" cells
%   [approach]_basin_certain_mean     - mean of the basin sizes, basins
%                                       formed by the "certain" cells
%   [approach]_basin_certain_median   - median of the basin sizes, basins
%                                       formed by the "certain" cells
%   [approach]_basin_certain_max      - maximum of the basin sizes, basins
%                                       formed by the "certain" cells
%   [approach]_basin_certain_std      - standard deviation of the basin 
%                                       sizes, basins formed by the 
%                                       "certain" cells
%   [approach]_basin_certain_sum      - sum of the basin sizes, basins
%                                       formed by the "certain" cells
%   [approach]_basin_uncertain_min    - minimum of the basin sizes, basins
%                                       formed by the "certain" and 
%                                       "uncertain" cells
%   [approach]_basin_uncertain_mean   - minimum of the basin sizes, basins
%                                       formed by the "certain" and 
%                                       "uncertain" cells
%   [approach]_basin_uncertain_median - median of the basin sizes, basins
%                                       formed by the "certain" and 
%                                       "uncertain" cells
%   [approach]_basin_uncertain_max    - maximum of the basin sizes, basins
%                                       formed by the "certain" and 
%                                       "uncertain" cells
%   [approach]_basin_uncertain_std    - standard deviation of the basin 
%                                       sizes, basins formed by the 
%                                       "certain" and "uncertain" cells
%   [approach]_basin_uncertain_sum    - sum of the basin sizes, basins
%                                       formed by the "certain" and 
%                                       "uncertain" cells
%   [approach]_best_attr_prob         - probability of finding attractor 
%                                       cell with the best objective value
%   [approach]_best_attr_no           - number of attractor cells with the 
%                                       best objective value

  if nargin < 3
    if nargin < 2
      help feature_gcm
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  lb = defopts(settings, 'lb', min(X));
  ub = defopts(settings, 'ub', max(X));
  blocks = defopts(settings, 'blocks', 2);
  metric = defopts(settings, 'distance', 'euclidean');
  dist_param = defopts(settings, 'dist_param', defMetricParam(metric, X));
  approach = defopts(settings, 'approach', {'min', 'mean', 'near'});
  min_fun = defopts(settings, 'minimize', true);
  
  nApproaches = numel(approach);
  % minimization
  if ~min_fun
    y = -y;
  end
  
  % create grid of cells
  cmg = CMGrid(X, y, lb, ub, blocks);
  nCells = cmg.nCells;
  
  sumCells = prod(cmg.blocks);
  % warn in case of empty cells
  if nCells < sumCells
    warning(['%d out of %d cells (%0.2f%%) is empty.', ...
             'This may make generalized cell mapping features useless.'], ...
            sumCells - nCells, sumCells, (sumCells - nCells)/sumCells * 100)
  end
  
  % get cell-function values
  for a = 1:nApproaches
    % gain y-values according to approach
    switch approach{a}
      case 'min'
        % get cell minimums
        [~, y_min] = cmg.getCellMin();
        y = NaN(cmg.blocks);
        y(coor2ind(cmg.cellId, cmg.blocks)) = y_min;
      case 'mean'
        % get cell means
        y_mean = cmg.getCellMean();
        y = NaN(cmg.blocks);
        y(coor2ind(cmg.cellId, cmg.blocks)) = y_mean;
      case 'near'
        % get points nearest to cell centers
        y = cmg.getNearCtrGridPointY(metric, dist_param);
      otherwise
        error('Approach %s is not implemented', approach{a})
    end
    
    RQ = [];
    periodicCells = false(1, sumCells);
    for c = 1:sumCells
      cellId = ind2coor(c, cmg.blocks);
      % get neighbor sets according to f(c_i) <= f(c_j)
      [bc, pg] = getNeighbors(cellId, cmg.blocks, y);
      % compute probability to go from c_i to c_j
      p = zeros(1, sumCells);
      p(c) = 1;
      if isempty(bc)
        p(pg) = 1/numel(pg);
      else
        p(bc) = (y(c) - y(bc)) / (sum(y(c) - y(bc)));
      end
      % is cell periodic?
      if sum(p) == 1
        periodicCells(c) = true;
      else
        RQ(end+1, :) = p;
      end
    end
    
    % split RQ to R and Q
    % R = RQ(:, periodicCells);
    Q = RQ(:, ~periodicCells);
    % compute fundamental matrix
    N = inv(eye(size(Q)) - Q);
    % if inversion is not possible, calculate less precise equivalent 
    % matrix
    % if any(any(isinf(N)))
    %   N = eye(size(Q));
    %   precision = 50;
    %   for k = 1:precision
    %     N = N + Q^k;
    %   end
    % end
    % compute absorbing probability
    % B = N*R;
    
    % calculate features
    ft.([approach{a}, '_attractors']) = size(N, 2);
    ft.([approach{a}, '_pcells']) = sum(periodicCells) / sumCells;
    ft.([approach{a}, '_tcells']) = sum(~periodicCells) / sumCells;
    ft.([approach{a}, '_uncertain']) = sum(sum(N ~= 0, 2) > 1) / sumCells;
    
    % compute probability of each basin of attraction
    basin_prob = sum(N) / sumCells;
    ft.([approach{a}, '_basin_prob_min']) = min(basin_prob);
    ft.([approach{a}, '_basin_prob_mean']) = mean(basin_prob);
    ft.([approach{a}, '_basin_prob_median']) = median(basin_prob);
    ft.([approach{a}, '_basin_prob_max']) = max(basin_prob);
    ft.([approach{a}, '_basin_prob_std']) = std(basin_prob);
    
    % compute basins of attraction
    basin_uncertain = calcBasinSize(N) / sumCells;
    basin_certain = calcBasinSize(N(sum(N ~= 0, 2) == 1, :)) / sumCells;
    ft.([approach{a}, '_basin_certain_min']) = min(basin_certain);
    ft.([approach{a}, '_basin_certain_mean']) = mean(basin_certain);
    ft.([approach{a}, '_basin_certain_median']) = median(basin_certain);
    ft.([approach{a}, '_basin_certain_max']) = max(basin_certain);
    ft.([approach{a}, '_basin_certain_std']) = std(basin_certain);
    ft.([approach{a}, '_basin_certain_sum']) = sum(basin_certain);
    
    ft.([approach{a}, '_basin_uncertain_min']) = min(basin_uncertain);
    ft.([approach{a}, '_basin_uncertain_mean']) = mean(basin_uncertain);
    ft.([approach{a}, '_basin_uncertain_median']) = median(basin_uncertain);
    ft.([approach{a}, '_basin_uncertain_max']) = max(basin_uncertain);
    ft.([approach{a}, '_basin_uncertain_std']) = std(basin_uncertain);
    ft.([approach{a}, '_basin_uncertain_sum']) = sum(basin_uncertain);
    
    % probability to find the best cell
    y_attr = y(~periodicCells);
    ft.([approach{a}, '_best_attr_prob']) = sum(basin_prob(y_attr == min(y_attr)));
    ft.([approach{a}, '_best_attr_no']) = sum(y_attr == min(y_attr)) / sumCells;
  end
  
  % ensure features to be non-empty in case of empty input
  if isempty(X) || isempty(y)
    ft = repStructVal(ft, @isempty, NaN, 'test');
  end
  
end

function [bc, pg] = getNeighbors(cellId, blocks, y)
% get bc and pg sets according to cellId and y-values

  bc = [];
  pg = [];
  cellInd = coor2ind(cellId, blocks);
  % find neighbors
  for d = 1 : numel(blocks)
    % left neighbor
    neighborId = cellId;
    neighborId(d) = neighborId(d) - 1;
    if neighborId(d) > 0 && neighborId(d) < blocks(d) + 1
      neighborInd = coor2ind(neighborId, blocks);
      if y(neighborInd) < y(cellInd)
        bc(end + 1) = coor2ind(neighborId, blocks);
      elseif y(neighborInd) == y(cellInd)
        pg(end + 1) = coor2ind(neighborId, blocks);
      end
    end
    % right neighbor
    neighborId = cellId;
    neighborId(d) = neighborId(d) + 1;
    if neighborId(d) > 0 && neighborId(d) < blocks(d) + 1
      neighborInd = coor2ind(neighborId, blocks);
      if y(neighborInd) < y(cellInd)
        bc(end + 1) = coor2ind(neighborId, blocks);
      elseif y(neighborInd) == y(cellInd)
        pg(end + 1) = coor2ind(neighborId, blocks);
      end
    end
  end
end

function basins = calcBasinSize(N)
% calculate basin size
  if isempty(N)
    basins = 0;
  else
    [~, max_id] = max(N);
    for i = 1:numel(max_id)
      basins(i) = sum(i == max_id);
    end
  end
end