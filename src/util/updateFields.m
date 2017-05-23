function s = updateFields(former, update, i, j)
  % Updates fields of the struct 'former' with the contents
  % of the respective fields from 'update'.
  % If the contents of some field in 'former' is cell array
  % or a numeric array, the new content from 'update' is written
  % at indices {i,j} or (i,j), i.e.
  %
  %   former.(fieldname){i,j} = update.(fieldname)
  %
  %   former.(fieldname)(i,j) = update.(fieldname)
  former_fields = fields(former)
  update_fields = fields(update)

  % cycle through all fields
  for f = 1:length(former_fields)
    fname = former_fields{i};
    % update only if there is such field in the 'former'...
    if (ismember(fname, update_fields))
      if (iscell(former.(fname)) && ~isempty(i) && ~isempty(j))
        % field is a cell array
        former.(fname){i,j} = update.(fname);
      elseif (isnumeric(former.(fname)) && (size(former.(fname),1) > 1 ...
          || size(former.(fname),2) > 1) && ~isempty(i) && ~isempty(j))
        % field is a non-trivial numeric array
        former.(fname)(i,j) = update.(fname);
      else
        former.(fname) = update.(fname);
      end
    end
  end
end
