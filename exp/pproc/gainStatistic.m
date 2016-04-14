function stats = gainStatistic(data, dimId, funcId, nInstances, averageDims, statistic)
% Returns cell array of means accross chosen dimensions for each function
%
% Input:
%   statistic - handle to statistic function | @mean, @median

  if nargin < 4 || isempty(nInstances)
    nInstances = 15;
    if nargin < 5
      averageDims = true;
      if nargin < 6
        statistic = @mean;
      end
    end
  end

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
        actualData = data{funcId(f), dimId(d)};
        useInstances = min([nInstances, size(actualData, 2)]);
        stats{f, d} = statistic(actualData(:, 1:useInstances), 2);
      end
    end
  end

end