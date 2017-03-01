function stats = gainStatistic(data, dimId, funcId, varargin)
% stats = gainStatistic(data, dimId, funcId, varargin) returns cell array 
% of means accross chosen dimensions for each function.
%
% Input:
%   data     - cell-array of data
%   dimId    - identifiers of dimensions | integer
%   funcId   - identifiers of functions | integer
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%
%     'Statistic'    - handle to statistic function | @mean, @median
%     'AverageDims'  - average accross dimensions | boolean
%     'MaxInstances' - maximum number of instances to use | integer
%     'SuppWarning'  - suppress warning if data in one function and 
%                      dimension are missing
%
% Output:
%   stats - means accross chosen dimensions for each function | cell-array

  stats = {};
  if nargin < 1
    help gainStatistic
    return
  end
  
  if isempty(data)
    warning('Data is empty')
    return
  end
  
  settings = settings2struct(varargin{:});
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
        if funcId(f)<= size(data, 1) && dimId(d) <= size(data,2) && ~isempty(data{funcId(f), dimId(d)})
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