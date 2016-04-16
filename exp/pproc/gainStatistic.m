function stats = gainStatistic(data, dimId, funcId, varargin)
% Returns cell array of means accross chosen dimensions for each function
%
% Input:
%   data
%   dimId
%   funcId
%   settings:
%     'Statistic' - handle to statistic function | @mean, @median
%     'AverageDims'
%     'MaxInstances'
%     'SuppWarning' - suppress warning if data in one function and 
%                     dimension are missing

  stats = {};
  if nargin < 1
    help gainStatistic
    return
  end
  
  if isempty(data)
    warning('Data is empty')
    return
  end
  
  if isstruct(varargin)
    settings = varargin;
  else
    % keep cells as cells due to struct command
    vararCellId = cellfun(@iscell, varargin);
    varargin(vararCellId) = {varargin(vararCellId)};
    settings = struct(varargin{:});
  end
  statistic = defopts(settings, 'Statistic', @mean);
  averageDims = defopts(settings, 'AverageDims', true);
  nInstances = defopts(settings, 'MaxInstances', 15);
  suppWarning = defopts(settings, 'SuppWarning', false);

  % cat dimensions if necessary
  dims = length(dimId);
  funcs = length(funcId);
  stats = cell(funcs, 1);
  if averageDims
    for f = 1:funcs
      funcData = [];
      for d = 1:dims
        actualData = data{funcId(f), dimId(d)};
        useInstances = min([nInstances, size(actualData, 2)]);
        funcData = [funcData, actualData(:, 1:useInstances)];
      end
      stats{f} = statistic(funcData, 2);
    end
  else
    for f = 1:funcs
      for d = 1:dims
        if ~isempty(data{funcId(f), dimId(d)})
          actualData = data{funcId(f), dimId(d)};
          useInstances = min([nInstances, size(actualData, 2)]);
          stats{f, d} = statistic(actualData(:, 1:useInstances), 2);
        else
          if ~suppWarning
            warning('Data in function %d and dimension %d is empty.', f, d)
          end
          stats{f, d} = [];
        end
      end
    end
  end

end