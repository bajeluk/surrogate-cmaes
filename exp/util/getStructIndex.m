function index = getStructIndex(origStruct, searchStruct)
% getStructIndex(origStruct, searchStruct) returns indices of searched 
% substructure in cell array of structures.
%
% Input:
%   origStruct   - cell array of structures for search
%   searchStruct - substructure to be searched in 'origStruct'
%
% Output:
%   index - vector of indices of searched substructure in 'origStruct'

  index = [];
  if nargin < 2
    help getStructIndex
    return
  end
  searchedFields = findFields(searchStruct);
  nFields = length(searchedFields);
  searchedValues = cell(1, nFields);
  % find all field values
  for j = 1:nFields 
    searchedValues{j} = eval(['searchStruct.', searchedFields{j}]);
  end
  for i = 1:length(origStruct)
    if ~isempty(origStruct{i})
      correctFields = true;
      % compare all needed fields
      for j = 1:nFields 
        try
          correctFields = correctFields && ...
                          all(isequal(eval(['origStruct{i}.', searchedFields{j}]), searchedValues{j}));
        catch err
          if strcmp(err.identifier, 'MATLAB:nonExistentField')
            % warning('Reference to non-existent field ''%s''.', searchedFields{j})
            % fprintf(2, 'Reference to non-existent field ''%s''.\n', searchedFields{j});
            correctFields = false;
          else
            throw(err)
          end
        end
      end
      if correctFields
        index(end+1) = i;
      end
    end
  end
end

function resFields = findFields(str)
% find recursive fields in structure

  resFields = {};
  strFields = fieldnames(str);
  nFields = length(strFields);
  for i = 1:nFields
    if isstruct(str.(strFields{i}))
      actualFields = findFields(str.(strFields{i}));
      resFields = [resFields, cellfun(@(x) [strFields{i}, '.', x], actualFields, 'UniformOutput', false)];
    else
      resFields = [resFields, strFields{i}];
    end
  end
end