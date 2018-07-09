function str = struct2str( S )
%STRUCT2STR Format a structure into a string.
  
  fs = cell(size(fields(S)));
  
  if ~numel(fs)
    str = '';
    return;
  end

  i = 1;
  for f = fieldnames(S)'
    fs{i} = sprintf('%s="%s"', f{:}, S.(f{:}));
    i = i + 1;
  end
  str = strjoin(fs, ', ');
end

