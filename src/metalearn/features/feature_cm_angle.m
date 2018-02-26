function ft = feature_cm_angle(X, y, settings)
% ft = FEATURE_CM_ANGLE(X, y, settings) returns cell mapping angle features 
% for dataset [X, y] according to settings.
%
% The cell mapping features discretizes the continuous input space
% utilizing a pre-defined number of blocks (cells) per dimension. Angle
% features extract information based on the location of the best and worst
% point of the cell considering the cell center.
%
% settings:
%   blocks   - number of cell blocks per dimension
%   lb       - lower bounds of the input space
%   minimize - binary flag stating whether the objective function should 
%   ub       - upper bounds of the input space
%
% Features:
%   dist_ctr2best_mean  - mean of distances from cell center to the best 
%                         point within the cell
%   dist_ctr2best_std   - standard deviation of distances from cell center 
%                         to the best point within the cell
%   dist_ctr2worst_mean - mean of distances from cell center to the worst 
%                         point within the cell
%   dist_ctr2worst_std  - standard deviation of distances from cell center 
%                         to the worst point within the cell
%   angle_mean          - mean of angles between the best, the worst, and
%                         the cell-center point within the cell
%   angle_std           - standard deviation of angles between the best, 
%                         the worst, and the cell-center point within the 
%                         cell
%   y_best2worst_mean   - mean of differences between the best and the 
%                         worst objective value per cell normalized by the
%                         (y_max - y_min) value
%   y_best2worst_std    - standard deviation of differences between the 
%                         best and the worst objective value per cell 
%                         normalized by the (y_max - y_min) value

  if nargin < 3
    if nargin < 2
      help feature_cm_angle
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end

  % parse settings
  min_fun = defopts(settings, 'minimize', true);
  lb = defopts(settings, 'lb', min(X));
  ub = defopts(settings, 'ub', max(X));
  blocks = defopts(settings, 'blocks', 1);
  
  if ~min_fun
    y = -y;
  end
  
  % create grid of cells
  cmCells = createCMGrid(X, y, lb, ub, blocks);
  nCells = numel(cmCells);
  
  % distances max - center in all cells
  distCentMax = cell2mat(arrayfun(@(i) cmCells(i).getDistCtr2Max, 1:nCells, ...
                                 'UniformOutput', false));
  % distances min - center in all cells
  distCentMin = cell2mat(arrayfun(@(i) cmCells(i).getDistCtr2Min, 1:nCells, ...
                                 'UniformOutput', false));
  % angle between min - center - max in all cells
  maxMinAngle = cell2mat(arrayfun(@(i) cmCells(i).getMaxMinAngle, 1:nCells, ...
                                 'UniformOutput', false)); 
  % difference y_max - y_min in all cells
  maxMinDiff = cell2mat(arrayfun(@(i) cmCells(i).getMaxMinDiff, 1:nCells, ...
                                 'UniformOutput', false));
  y_mnmx = minmax(y');
  
  % calculate features
  ft.dist_ctr2best_mean = mean(distCentMax);
  ft.dist_ctr2best_std  =  std(distCentMax);
  ft.dist_ctr2worst_mean = mean(distCentMin);
  ft.dist_ctr2worst_std  =  std(distCentMin);
  ft.angle_mean = mean(maxMinAngle);
  ft.angle_std  =  std(maxMinAngle);
  ft.y_best2worst_mean = mean(maxMinDiff/(y_mnmx(1)-y_mnmx(2)));
  ft.y_best2worst_std  =  std(maxMinDiff/(y_mnmx(1)-y_mnmx(2)));
end