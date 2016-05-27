function colors = getAlgColors(n)
% colors = getAlgColors(n) return colors for algorithms in report 
% publishing

%TODO: colors should be generated through some deterministic sequence
%      returning values far enough from each other (and from reference
%      algorithms colors)
  
  colors = [];
  if nargin < 1
    help getAlgColors
    return
  end
  
  hue_start = 3;
  hue = hue_start./(hue_start+1 : n+hue_start);
  hue = mod(cumsum(hue'), 1);
%   hue = rand(n, 1);
  
  sat_base = [1; 1; 0.5];
  val_base = [1; 0.5; 1];
  sat = repmat(sat_base, ceil(n/length(sat_base)), 1);
  val = repmat(val_base, ceil(n/length(val_base)), 1);
  sat = sat(1:n);
  val = val(1:n);
  colors = 255*hsv2rgb([hue, sat, val]);
end