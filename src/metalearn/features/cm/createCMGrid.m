function [cmCells, cmId] = createCMGrid(X, y, lb, ub, blocks)
% createCMGrid creates grid for calculating cell mapping metafeatures.
% 
% [cmCells, cmId] = createCMGrid(X, y, lb, ub, blocks)
%
% Input:
%   X - input features
%   y - input values
%   lb - lower bounds of input space
%   ub - upper bounds of input space
%   blocks - number of blocks (cells) per dimension 
%
% Output:
%   cmCells - array of mapping cells
%   cmId    - cell coordinates in discretized input space
%
% Example:
%   
%   X = rand(30, 3);
%   y = rand(30, 1);
%   lb = [0, 0, 0];
%   ub = [1, 1, 1];
%   blocks = [2, 3, 2];
%
%   [cmCells, cmId] = createCMGrid(X, y, lb, ub, blocks)
% 
%   cmCells = 
% 
%     1×12 CMCell array with properties:
% 
%       center
%       dim
%       lb
%       ub
%       minx
%       maxx
%       miny
%       maxy
%       X
%       y
% 
% 
%   cmId =
% 
%        1     2     1     2     1     2     1     2     1     2     1     2
%        1     1     2     2     3     3     1     1     2     2     3     3
%        1     1     1     1     1     1     2     2     2     2     2     2

  if nargin < 5
    if nargin < 4
      if nargin < 3
        if nargin < 2
          help createCMGrid
          if nargout > 0
            cmCells = CMCell();
            cmId = [];
          end
          return
        end
        lb = min(X);
      end
      ub = max(X);
    end
    blocks = 1;
  end
  
  [nData, dim] = size(X);
  
  % checkout blocks settings input
  lb = checkBlockVal(lb, 'lb', dim);
  assert(all(lb <= min(X)), 'Some points are out of lower bounds. Check your settings.')
  ub = checkBlockVal(ub, 'ub', dim);
  assert(all(ub >= max(X)), 'Some points are out of upper bounds. Check your settings.')
  blocks = checkBlockVal(blocks, 'blocks', dim);
  assert(all(blocks > 0 & mod(blocks, 1) == 0), 'Block numbers has to be natural.')
  
  % init
  blockLB = cell(1, dim);
  blockUB = cell(1, dim);
  pointCellId = zeros(nData, dim);
  cellIdVec = cell(1, dim);
  
  % calculate centers
  blockSize = (ub-lb)./blocks;
  for d = 1 : dim
    % lower bounds
    blockLB{d} = lb(d) + (0:(blocks(d)-1)) * blockSize(d);
    % upper bounds
    blockUB{d} = blockLB{d} + blockSize(d);
    % for each point find containing cells 
    pointCellId(:, d) = arrayfun(@(i) find(X(i, d) <= blockUB{d}, 1, 'first'), 1:nData);
    % create id vectors for cells
    cellIdVec{d} = 1:blocks(d);
  end
  
  % all possible cell coordinate combinations
  cmId = combvec(cellIdVec{:});
  
  % create cells
  for c = 1:prod(blocks)
    % find points related to actual cell
    actualPointId = all((repmat(cmId(:, c)', nData, 1) == pointCellId), 2);
    cellX = X(actualPointId, :);
    celly = y(actualPointId);
    cellSet.dim = dim;
    cellSet.lb = arrayfun(@(d) blockLB{d}(cmId(d, c)), 1:dim);
    cellSet.ub = arrayfun(@(d) blockUB{d}(cmId(d, c)), 1:dim);
    cmCells(c) = CMCell(cellX, celly, cellSet);
  end

end

function val = checkBlockVal(val, name, dim)
% check value of block setting
  if numel(val) == 1
    val = ones(1, dim) * val;
  elseif numel(val) ~= dim
    error('%s length differs from dimension', name)
  end
end