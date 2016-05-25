function colors = getAlgColors(n)
% colors = getAlgColors(n) return colors for algorithms in report 
% publishing

%TODO: colors should be generated through some deterministic sequence
%      returning values far enough from each other (and from reference
%      algorithms colors)
  
  hue = rand(n, 1);
  colors = hsv2rgb([hue, ones(n,1), ones(n,1)]);
end