function means = gainMeans(data, dimId, funcId, nInstances)
% Returns cell array of means accross chosen dimensions for each function

  if nargin < 4
    nInstances = 15;
  end

  % cat dimensions if necessary
  Dims = length(dimId);
  funcs = length(funcId);
  means = cell(funcs,1);
  if Dims > 1
    for f = 1:funcs
      funcData = [];
      for d = 1:Dims
        actualData = data{funcId(f),dimId(d)};
        useInstances = min([nInstances,size(actualData,2)]);
        funcData = [funcData,actualData(:,1:useInstances)];
      end
      means{f} = mean(funcData,2);
    end
  else
    for f = 1:funcs
      actualData = data{funcId(f),dimId};
      useInstances = min([nInstances,size(actualData,2)]);
      means{f} = mean(actualData(:,1:useInstances),2);
    end
  end

end