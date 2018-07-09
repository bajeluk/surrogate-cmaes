function node2dot(fid, i, lbl, props)
%NODE2DOT Print a node of a classification tree in dot format.
  lbl = ['label="' lbl, '"'];
  props = struct2str(props);
  parts = {lbl, props};
  parts = parts(cellfun(@(a) ~isempty(a), parts));
  fprintf(fid, '\t%d[%s]\n', i, strjoin(parts, ', '));
end

