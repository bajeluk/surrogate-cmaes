function ft = feature_basic(X, y, dataSettings)
% ft = FEATURE_BASIC(X, y, dataSettings) returns basic features for dataset
% [X, y]. Features provide obvious information of the initial design.
%
% dataSettings:
%   lb     - lower bounds
%   ub     - upper bounds
%   blocks - numbers of blocks per dimension
%   filled - indicator of filled blocks
%
% Features:
%   dim           - dimension of dataset
%   observations  - number of observations
%   lower_min     - minimum of lower bounds
%   lower_max     - maximum of lower bounds
%   upper_min     - minimum of upper bounds
%   upper_max     - maximum of upper bounds
%   objective_min - minimum of y values
%   objective_max - maximum of y values
%   blocks_min    - minimal number of blocks per dimension
%   blocks_max    - maximal number of blocks per dimension
%   cells_total   - total number of cells
%   cells_fille   - number of filled cells
%   minimize_fun  - binary flag stating whether the objective function
%                   should be minimized

  if nargin < 3
    if nargin < 2
      help feature_basic
      if nargout > 0
        ft = struct();
      end
      return
    end
    % TODO: proper settings
    dataSettings.lb = min(X);
    dataSettings.ub = max(X);
    dataSettings.blocks = 1;
    dataSettings.filled = ~isempty(X);
  end
  
  ft.dim = size(X, 2);
  ft.observations = numel(y);
  ft.lower_min = min(dataSettings.lb);
  ft.lower_max = max(dataSettings.lb);
  ft.upper_min = min(dataSettings.ub);
  ft.upper_max = max(dataSettings.ub);
  ft.objective_min = min(y);
  ft.objective_max = max(y);
  % TODO: write following lines properly
  ft.blocks_min = min(dataSettings.blocks);
  ft.blocks_max = max(dataSettings.blocks);
  ft.cells_total = prod(dataSettings.blocks);
  ft.cells_filled = sum(dataSettings.filled);
  ft.minimize_fun = false;

end