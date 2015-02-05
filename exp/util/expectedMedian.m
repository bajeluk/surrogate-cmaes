function m = expectedMedian(x, varargin)

  if ((nargin > 1) && varargin{1} == 2)
    % along dimension 2
    x = x';
    m = nan(size(x,2), 1);
  else
    % along dimension 1
    m = nan(1, size(x,2));
  end
  % penalize the NaN's? Default is yes.
  penalize = ~((nargin > 2) && ~varargin{2});

  for i = 1:size(x,2)
    nans = isnan(x(:,i));
    if (sum(~nans) > 0)
      % calculate median for non-NaN's
      m(i) = median(x(~nans,i));
      % divide the number by the ratio of NaN's
      if (sum(nans) > 0  &&  penalize)
        m(i) = m(i) * length(nans) / sum(~nans);
      end
    end
  end
end
