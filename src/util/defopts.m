function value = defopts(structOpts, fieldName, defValue)
  % return either value specified in structOpts.(fieldName)  (if specified)
  % or  'defValue'  (otherwise)
  if (isfield(structOpts, fieldName) && ~isempty(structOpts.(fieldName)))
    value = structOpts.(fieldName);
  else
    value = defValue;
  end
end
