function props = meta_leaf_fmt(cls)
  if strcmp(cls, 'g')
    props = struct('shape', 'box', 'fillcolor', 'Green', 'style', 'filled');
  elseif strcmp(cls, 'b')
    props = struct('shape', 'box', 'fillcolor', 'LightBlue', 'style', 'filled');
  else
    props = struct('shape', 'box');
  end
end

