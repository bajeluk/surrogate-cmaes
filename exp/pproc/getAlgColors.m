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
    ];
  
  def_alg_color_base = [...
     22,  22, 138; ... % dark blue    | CMA-ES
    255, 225,   0; ... % yellow       | lmm-CMA-ES
    100, 149, 237; ... % light blue   | BIPOP-saACMES
    178,  34,  34; ... % bloody red   | S-CMA-ES
    154, 205,  50; ... % some green   | DTS-CMA-ES
    255, 155,   0; ... % light orange | SMAC
    ];
  
  def_alg_names = {'cmaes', 'lmmcmaes', 'saacmes', 'scmaes', 'dtscmaes', 'smac'};
 
  if ischar(colId)
    colId = {colId};
  end
  
  % return default algorithm colors
  if iscell(colId)
    colors = zeros(length(colId), 3);
    
    algId = cellfun(@(x) find(strcmp(x, def_alg_names)), colId, 'UniformOutput', false);
    noAlg = cellfun(@isempty, algId);
    if any(noAlg)
      warning('Algorithms %s have no default color. Returned colors are from the predefined spectrum.', ...
              strjoin(colId(noAlg), ', '));
      colors_rep = repmat(color_base, ceil(sum(noAlg)/length(color_base)), 1);
      colors( noAlg, :) = colors_rep(1:sum(noAlg), :);
    end
    colors(~noAlg, :) = def_alg_color_base(cell2mat(algId), :);
    
  % return common colors
  else
    max_color = max(colId);
    colors = repmat(color_base, ceil(max_color/length(color_base)), 1);
    colors = colors(colId, :);
  end
end