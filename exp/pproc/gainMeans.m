function means = gainStatistic(data, dimId, funcId, nInstances, averageDims, statistics)
% Returns cell array of means accross chosen dimensions for each function

  if nargin < 4
    nInstances = 15;
    if nargin < 5
      averageDims = true;
    end
  end

  % cat dimensions if necessary
  dims = length(dimId);
  funcs = length(funcId);
  means = cell(funcs, 1);
  if averageDims
    for f = 1:funcs
      funcData = [];
      for d = 1:dims
        actualData = data{funcId(f), dimId(d)};
        useInstances = min([nInstances, size(actualData, 2)]);
        funcData = [funcData, actualData(:, 1:useInstances)];
      end
      means{f} = mean(funcData, 2);
    end
  else
    for f = 1:funcs
      for d = 1:dims
        actualData = data{funcId(f), dimId(d)};
        useInstances = min([nInstances, size(actualData, 2)]);
        means{f, d} = mean(actualData(:, 1:useInstances), 2);
      end
    end
  end

end