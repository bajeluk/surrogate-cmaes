function y = mapColumns(x, fun, varargin)
% apply function to columns (or rows if 3rd argument == 2),
% but ignore NaNs

  if ((nargin > 2) && varargin{1} == 2)
    % along dimension 2
    x = x';
    y = nan(size(x,2), 1);
  else
    % along dimension 1
    y = nan(1, size(x,2));
  end
  
  for i = 1:size(x,2)
    nans = isnan(x(:,i));
    if (sum(~nans) > 0)
      y(i) = fun(x(~nans,i));
    end
  end
end
