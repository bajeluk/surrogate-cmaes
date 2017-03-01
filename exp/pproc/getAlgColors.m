function colors = getAlgColors(colId)
% colors = getAlgColors(colId) return colors with ID colID for algorithms 
% in report publishing.
%
% Example:
%   colors = getAlgColors(1:4)
%
%   colors =
%
%      255   225     0
%      175   248     6
%       12   240   248
%      248    86     6
  
  if nargout > 1
    colors = [];
  end
  if nargin < 1
    help getAlgColors
    return
  end

  color_base = [...
    175, 248,   6; ... % light green
     12, 240, 248; ... % azure (almost cyan)
    248,  86,   6; ... % dark orange
      0, 255, 120; ... % alien green
    154,  22, 106; ... % light violet
    255,   0,   0; ... % red
     20, 205,  16; ... % green
     97, 143, 163; ... % metal grey
     88,  51, 138; ... % dark violet
     13, 129,  10; ... % dark green
    190, 183,  58; ... % no idea green
    248,   6, 246; ... % pink-violet
     27, 147, 134; ... % kerosene (blue-green)
      0,   0,   0; ... % black
    255, 112, 168; ... % pink (chewing gum)
    149, 149, 149; ... % light grey
    149, 128,  78; ... % khaki
    163,  97,  97; ... % almost brown
     77,  77,  77 ...  % dark grey
    %  36, 140, 248; ... % light blue   | BIPOP-saACMES
    %  22,  22, 138; ... % dark blue    | CMA-ES
    % 154, 205,  50; ... % some green   | DTS-CMA-ES
    % 178,  34,  34; ... % bloody red   | S-CMA-ES
    % 255, 155,   0; ... % light orange | SMAC
    % 255, 225,   0; ... % yellow       | lmm-CMA-ES
    ];
  
  max_color = max(colId);
  colors = repmat(color_base, ceil(max_color/length(color_base)), 1);
  colors = colors(colId, :);
end