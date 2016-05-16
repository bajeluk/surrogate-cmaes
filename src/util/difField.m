function [df, dVals] = difField(structCell)
% return different field values and its fieldnames

  df = {};
  dVals = {};
  refStruct = structCell{1};
  for s = 1:length(structCell)
    
    %TODO: use function subfields from printStructure
    sNames = fieldnames(structCell{s});
    for n = 1:length(sNames)
      if ~isfield(refStruct, sNames{n}) || ~isequal(refStruct.(sNames{n}), structCell{s}.(sNames{n}))
        df{end+1} = sNames{n};
        dVals{end+1} = structCell{s}.(sNames{n});
      end
    end
  end
end