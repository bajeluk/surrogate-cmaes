function ft = feature_basic(X, y, settings)
% ft = FEATURE_BASIC(X, y, settings) returns basic features for dataset
% [X, y] according to settings. Features provide obvious information of the
% initial design.
%
% settings:
%   lb       - lower bounds | default: min(X)
%   ub       - upper bounds | default: max(X)
%   blocks   - numbers of blocks per dimension | default: 1
%   filled   - indicator of filled blocks | default: true
%   minimize - binary flag stating whether the objective function should be
%              minimized | default: true
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
    settings = struct();
  end
  
  % TODO: proper settings
  lb = defopts(settings, 'lb', min(X));
  ub = defopts(settings, 'ub', max(X));
  blocks = defopts(settings, 'blocks', 1);
  filled = defopts(settings, 'filled', ~isempty(X));
  min_fun = defopts(settings, 'minimize', true);
  
  ft.dim = size(X, 2);
  ft.observations = numel(y);
  ft.lower_min = min(lb);
  ft.lower_max = max(lb);
  ft.upper_min = min(ub);
  ft.upper_max = max(ub);
  ft.objective_min = min(y);
  ft.objective_max = max(y);
  % TODO: write following lines properly
  ft.blocks_min = min(blocks);
  ft.blocks_max = max(blocks);
  ft.cells_total = prod(blocks);
  ft.cells_filled = sum(filled);
  ft.minimize_fun = min_fun;

end