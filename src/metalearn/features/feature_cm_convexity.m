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
% (soft/hard) convexity or concavity. (Kerschke et al., 2017)
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
  % default dist_param value
  switch metric
    case 'mahalanobis'
      def_param = nancov(X);
    case 'minkowski'
      def_param = 2;
    case 'seuclidean'
      def_param = nanstd(X);
    otherwise
      def_param = [];
  end
  dist_param = defopts(settings, 'dist_param', def_param);
  
  % checkout number of cells per dimension
  if any(blocks < 3)
    warning('The minimal number of cells per dimension is 3.')
    ft = struct();
  end
  
  % create grid of cells
  [cmCells, cmId] = createCMGrid(X, y, lb, ub, blocks);
  nCells = numel(cmCells);
  dim = cmCells(1).dim;
  
  % calculate points nearest to cell center in all cells
  y = NaN(nCells, 1);
  for i = 1:nCells
    [~, y_act] = cmCells(i).getNearCtrPoint(metric, dist_param);
    if isempty(y_act)
      y(i) = NaN;
    else
      y(i) = y_act;
    end
  end
  
  % warn in case of empty cells
  if any(isnan(y))
    warning('%d out of %d cells (%0.2f%%) is empty. This may affect the results.', ...
            sum(isnan(y)), nCells, sum(isnan(y))/nCells * 100)
  end
  
  % calculate convexity and concavity per dimension
  conc_soft = 0;
  conv_soft = 0;
  conc_hard = 0;
  conv_hard = 0;
  nComparisons = 0;
  for d = 1 : dim
    for b = 2 : blocks(d) - 1
      % identify neighbor values
      left_y   = y(cmId(d, :) == b - 1);
      center_y = y(cmId(d, :) == b);
      right_y  = y(cmId(d, :) == b + 1);
      % compare for soft/hard convexity/concavity
      % ids are periodical => not necessary to checkout other positions
      conc_soft = conc_soft + sum(center_y > (left_y + right_y)/2);
      conv_soft = conv_soft + sum(center_y < (left_y + right_y)/2);
      conc_hard = conc_hard + sum(center_y > max(left_y, right_y));
      conv_hard = conv_hard + sum(center_y < min(left_y, right_y));
      % increase number of comparisons
      nComparisons = nComparisons + numel(center_y);
    end
  end
  
  % calculate features
  ft.concave_soft = conc_soft/nComparisons;
  ft.concave_hard = conc_hard/nComparisons;
  ft.convex_soft  = conv_soft/nComparisons;
  ft.convex_hard  = conv_hard/nComparisons;
end