function stars = decorateKendall(kendall)
  % transfer kendall rho coefficient into a string
  %
  % [ ***** ] ... [ *     ] [   o   ] [ -     ]  ... [ ----- ]
  %    1.0    ... 0.01-0.19    0.0    -0.01-0.19 ...   -1.0
  kendallInStars = floor(abs((kendall) * 5));
  if (kendallInStars == 0)
    stars = '[   o   ]';
  elseif (isnan(kendallInStars))
    stars = '[! NaN !]';
  else
    if (kendall > 0)
      mark = '*';
    else
      mark = '-';
    end
    space = ' ';
    stars = sprintf('[ %s%s ]', mark(ones(1,kendallInStars)), space(ones(1,5-kendallInStars)));
  end
end
