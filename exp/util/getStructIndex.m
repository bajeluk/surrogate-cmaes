function index = getStructIndex(origStruct, searchStruct)
% Returns indices of searched substructure in cell array of structures.
%
% Input:
%   origStruct   - cell array of structures for search
%   searchStruct - substructure to be searched in 'origStruct'
%
% Output:
%   index - vector of indices of searched substructure in 'origStruct'

  index = [];
  searchedFields = fieldnames(searchStruct);
  nFields = length(searchedFields);
  searchedValues = cell(1, nFields);
  for j = 1:nFields % find all field values
    searchedValues{j} = getfield(searchStruct, searchedFields{j});
  end
  for i = 1:length(origStruct)
    correctFields = true;
    for j = 1:nFields % compare all needed fields
      correctFields = correctFields && all(getfield(origStruct{i}, searchedFields{j}) == searchedValues{j});
    end
    if correctFields
      index(end+1) = i;
    end
  end
end