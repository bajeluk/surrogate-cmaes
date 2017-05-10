function settings = settings2struct(varargin)
% settings = settings2struct(settings) parses settings from cell-array in
% name - value pairs to structure format.
%
% Example:
%   All the following commands give the same result:
%
%   settings2struct('a', 1, 'b', 2)
%   settings2struct({'a', 1, 'b', 2})
%   settings2struct({'a', 1}, {'b', 2})
%   ans = 
%   
%       a: 1
%       b: 2

  if nargin < 1 || isempty(varargin) || isempty(varargin{1})
    settings = struct();
    return
  end

  % struct input is immediately returned
  if isstruct(varargin{1})
    settings = varargin{1};
    return
  elseif (iscell(varargin{1}) && isstruct(varargin{1}{1}))
    settings = varargin{1}{1};
    return
  end
  % multiple cell input
  if all(cellfun(@iscell, varargin))
    varargin = [varargin{:}];
  end
  % necessary name-value conditions
  assert(mod(length(varargin), 2) == 0, 'Number of cell array elements has to be even')
  assert(all(cellfun(@ischar, varargin(1:2:end-1))), 'All name parameters has to be string')

  % keep cells as cells due to struct command
  vararCellId = cellfun(@iscell, varargin);
  varargin(vararCellId) = cellfun(@(x) {x}, varargin(vararCellId), 'UniformOutput', false);
  settings = struct(varargin{:});
end