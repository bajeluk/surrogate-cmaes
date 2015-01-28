function [bbParams, sgParams, cmParams, nNonBbobValues] = getParamsFromIndex(id, bbAll, sgAll, cmAll)
% GETPARAMSFROMINDEX Generates struct arrays with parameters accrd. to exper. ID#
  bbNValues = cell2mat(structMap(bbAll, @(x) length(x.values)));
  sgNValues = cell2mat(structMap(sgAll, @(x) length(x.values)));
  cmNValues = cell2mat(structMap(cmAll, @(x) length(x.values)));

  nValues = [bbNValues sgNValues cmNValues];
  nParams = length(nValues);
  totalCombs = prod(nValues);

  assert(id <= totalCombs && id > 0, 'getParamsFromIndex(): the specified index is greater than total # combinations.');

  orders = size(nValues);
  for i = 1:nParams
    orders(i) = prod(nValues(i:end));
  end
  % orders =  [540   180    60    60    60    30    10    10    10     2];
  % nValues =   [3     3     1     1     2     3     1     1     5     2];
  orders = [orders 1];

  % The number of different parameter-values combinations except BBOB settings
  nNonBbobValues = prod(nValues((length(bbAll)+1):end));

  x = id-1;
  i = 1;

  % vector of integer-parameter-values:
  paramIDs = zeros(size(nValues));

  while (x >= 0 && i <= nParams)
    div = floor(x / orders(i+1));
    paramIDs(i) = div + 1;
    x = mod(x, orders(i+1));
    i = i + 1;
  end

  pid = 0;

  % BBOB parameters
  bbParams = struct();
  for i = 1:length(bbAll)
    bbParams.(bbAll(i).name) = bbAll(i).values{paramIDs(i)};
  end
  pid = pid + length(bbAll);

  % Surrogate parameters
  sgParams = struct();
  for i = 1:length(sgAll)
    fieldnames = strsplit(sgAll(i).name, '.');
    switch length(fieldnames)
    case 1
      sgParams.(sgAll(i).name) = sgAll(i).values{paramIDs(pid+i)};
    case 2
      sgParams.(fieldnames{1}).(fieldnames{2}) = sgAll(i).values{paramIDs(pid+i)};
    case 3
      sgParams.(fieldnames{1}).(fieldnames{2}).(fieldnames{3}) = sgAll(i).values{paramIDs(pid+i)};
    otherwise
      % join the names with underscores
      sgParams.(strjoin(C,'_')) = sgAll(i).values{paramIDs(pid+i)};
    end
  end
  pid = pid + length(sgAll);

  % CMA-ES parameters
  cmParams = struct();
  for i = 1:length(cmAll)
    cmParams.(cmAll(i).name) = cmAll(i).values{paramIDs(pid+i)};
  end
end
