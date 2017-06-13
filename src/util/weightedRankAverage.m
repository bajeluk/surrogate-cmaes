function awx = weightedRankAverage(x, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  x = x(:);
  if (nargin >= 2)
    if (isnumeric(varargin{1}))
      weights = varargin{1}(:);
    elseif (ischar(varargin{1}) && strcmpi(varargin{1}, 'reverse'))
      weights = (1:length(x))';
    end
  else
    weights = (length(x):-1:1)';
  end
  
  sx = sort(x);
  awx = sum(sx .* weights) / sum(weights);
end
