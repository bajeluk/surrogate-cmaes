function multiFieldNames = getFieldsWithMultiValues(modelOpts)
% Get names of fields with more than one value in its cell-array fieldvalue
  assert(isstruct(modelOpts));

  fnames = fieldnames(modelOpts);
  % multiFieldNames = cell(1, length(fnames));
  multiFieldNames = {};
  for i = 1:length(fnames)
    fname = fnames{i};
    if (iscell(modelOpts.(fname)) && length(modelOpts.(fname)) > 1)
      multiFieldNames{end+1} = fname;
    end
  end
end