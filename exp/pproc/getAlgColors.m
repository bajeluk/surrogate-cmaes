function colors = getAlgColors(n)
% colors = getAlgColors(n) return colors for algorithms in report 
% publishing

%TODO: colors should be generated through some deterministic sequence
%      returning values far enough from each other (and from reference
%      algorithms colors)
  colors = randi(256, n, 3) - 1;
end