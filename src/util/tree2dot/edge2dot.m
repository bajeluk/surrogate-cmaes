function edge2dot(fid, i, j, lbl, props)
%EDGE2DOT Print an edge of a classification tree in dot format.
  lbl = ['label="' lbl, '"'];
  props = struct2str(props);
  parts = {lbl, props};
  parts = parts(cellfun(@(a) ~isempty(a), parts));
  fprintf(fid, '\t%d -> %d[%s]\n', i, j, strjoin(parts, ', '));
end
