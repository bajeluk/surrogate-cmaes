function ft = feature_cm_gradhomo(X, y, settings)
% ft = FEATURE_CM_GRADHOMO(X, y, settings) returns cell mapping gradient
% homogeneity features for dataset [X, y] according to settings.
%
% The cell mapping features discretizes the continuous input space
% utilizing a pre-defined number of blocks (cells) per dimension. Gradient 
% homogeneity features aggregate the cell-wise information on the gradients
% between each point of a cell and its corresponding nearest neighbor. For 
% each point within a cell, we compute the gradient towards its nearest 
% neighbor, normalize it, point it towards the better one of the two 
% neighboring points and afterwards sum up all the normalized gradients 
% (per cell). (Kerschke et al., 2017)
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
%   grad_mean - mean of homogeneity gradients accross all cells
%   grad_std  - standard deviation of homogeneity gradients accross all 
%               cells

  if nargin < 3
    if nargin < 2
      help feature_cm_gradhomo
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
  blocks = defopts(settings, 'blocks', 1);
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
  
  % create grid of cells
  cmCells = createCMGrid(X, y, lb, ub, blocks);
  nCells = numel(cmCells);
  
  % calculate gradient homogeneity in all cells
  gradHomo = cell2mat(arrayfun(@(i) cmCells(i).getGradHomogeneity(metric, dist_param), ...
                               1:nCells, 'UniformOutput', false));
                             
  % checkout number of filled cells
  if nCells > numel(gradHomo)
    warning('%d out of %d cells (%0.2f%%) contain less than three observations.', ...
            nCells - numel(gradHomo), nCells, (nCells - numel(gradHomo))/nCells * 100)
  end
  
  % calculate features
  ft.grad_mean = mean(gradHomo);
  ft.grad_std  =  std(gradHomo);
end