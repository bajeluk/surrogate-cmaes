function ft = feature_cm_convexity(X, y, settings)
% ft = FEATURE_CM_CONVEXITY(X, y, settings) returns cell mapping convexity
% features for dataset [X, y] according to settings.
%
% The cell mapping features discretizes the continuous input space
% utilizing a pre-defined number of blocks (cells) per dimension. 
% Convexity features aggregate the (estimated) convexity based on 
% representative observations of successive cells. Each cell is 
% represented by the sample observation that is located closest to the 
% corresponding cell center. We then compare – for each combination of 
% three neighboring cells – their objective values: if the objective value 
% of the middle cell is below/above the mean of the objective values of the 
% outer two cells, we have an indication for (soft) convexity/concavity. 
% If the objective value of the middle cell even is the lowest/biggest 
% value of the three neighboring cells, we have an indication for (hard) 
% convexity/concavity. Averaging these values across all combinations of 
% three neighboring cells, results in the estimated “probabilities” for 
% (soft/hard) convexity or concavity. (Kerschke et al., 2014)
%
% settings:
%   blocks     - number of cell blocks per dimension
%   distance   - distance metric (similar to pdist function) | default:
%                'euclidean'
%   dist_param - additional parameter to distance (similar to pdist 
%                function)
%   lb         - lower bounds of the input space
%   ub         - upper bounds of the input space
%
% Features:
%   concave_soft - estimated probability of soft concavity
%   concave_hard - estimated probability of hard concavity
%   convex_soft  - estimated probability of soft convexity
%   convex_hard  - estimated probability of hard convexity

  if nargin < 3
    if nargin < 2
      help feature_cm_convexity
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
  blocks = defopts(settings, 'blocks', 3);
  metric = defopts(settings, 'distance', 'euclidean');
  dist_param = defopts(settings, 'dist_param', defMetricParam(metric, X));
  
  % checkout number of cells per dimension
  if any(blocks < 3)
    warning('The minimal number of cells per dimension is 3.')
    ft = struct();
  end
  
  % create grid of cells
%   [cmCells, cmId] = createCMGrid(X, y, lb, ub, blocks);
  cmg = CMGrid(X, y, lb, ub, blocks);
  nCells = cmg.nCells;
  sumCells = prod(blocks);
  dim = cmg.dim;
  
  % calculate objective values of points nearest to cell center in all 
  % cells
  [~, yc] = cmg.getNearCtrPoint(metric, dist_param);
  
  % warn in case of empty cells
  if nCells < sumCells
    warning('%d out of %d cells (%0.2f%%) is empty. This may affect the results.', ...
            sumCells - nCells, sumCells, (sumCells - nCells)/sumCells * 100)
  end
  
  % calculate convexity and concavity per dimension
  conc_soft = 0;
  conv_soft = 0;
  conc_hard = 0;
  conv_hard = 0;
  nComparisons = 0;
  for d = 1 : dim
    % find if any non-empty cells are neighbours in this dimension
    dimCellId = cmg.cellId;
    dimCellId(:, d) = [];
    [~, ~, ud] = unique(dimCellId, 'rows');
    
    % potential neighbours
    for n = 1:max(ud)
      % ids of cells potential neighbors
      oneDimCellId = find(n == ud);
      % at least three cells in one dimension having all other dimensions
      % the same
      if numel(oneDimCellId) > 2
        % coordinates in recent dimension
        oneDimCoor = cmg.cellId(oneDimCellId, d);
        % search center cells one by one
        for centerCellId = 2 : numel(oneDimCellId) - 1
          % three neighboring cells
          if oneDimCoor(centerCellId) - 1 == oneDimCoor(centerCellId - 1) && ...
             oneDimCoor(centerCellId) + 1 == oneDimCoor(centerCellId + 1)
            left_y   = yc(oneDimCellId(centerCellId - 1));
            center_y = yc(oneDimCellId(centerCellId));
            right_y  = yc(oneDimCellId(centerCellId + 1));
            % compare for soft/hard convexity/concavity
            % ids are periodical => not necessary to checkout other positions
            conc_soft = conc_soft + sum(center_y > (left_y + right_y)/2);
            conv_soft = conv_soft + sum(center_y < (left_y + right_y)/2);
            conc_hard = conc_hard + sum(center_y > max(left_y, right_y));
            conv_hard = conv_hard + sum(center_y < min(left_y, right_y));
          end
        end
      end
    end

    % increase the number of comparisons including comparisons not
    % performed due to the lack of non-empty cells
    nComparisons = nComparisons + (blocks(d) - 2) * prod(blocks((1:dim) ~= d));
  end
  
  % calculate features
  ft.concave_soft = conc_soft/nComparisons;
  ft.concave_hard = conc_hard/nComparisons;
  ft.convex_soft  = conv_soft/nComparisons;
  ft.convex_hard  = conv_hard/nComparisons;
end