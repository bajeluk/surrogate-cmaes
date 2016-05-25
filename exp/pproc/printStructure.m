function str = printStructure(structure, FID, varargin)
% printStructure(structure, FID, option) prints all fields and values of 
% 'structure' to file 'FID'
%
% Input:
%   structure - structure to print
%   FID       - file identifier or filename | integer or string
%   settings  - pairs of property (string) and value or struct with 
%               properties as fields:
%                 'StructName' - name of structure to be printed in result
%                 'Format'     - format of resulting string:
%                                  'value'  - returns only field values
%                                  'field'  - does not return structure 
%                                             name

%TODO: include FID to varargin:
%         printStructure(structure)
%         printStructure(structure, FID)
%         printStructure(structure, varargin)
%         printStructure(structure, FID, varargin)

  str = [];
  if nargin < 2
    FID = 1;
    if nargin < 1
      help printStructure
      return
    end
  end
  
  if isstruct(varargin)
    settings = varargin;
  else
    % keep cells as cells due to struct command
    vararCellId = cellfun(@iscell, varargin);
    varargin(vararCellId) = {varargin(vararCellId)};
    settings = struct(varargin{:});
  end
  
  % parse settings
  structureName = defopts(settings, 'StructName', inputname(1));
  if isempty(structureName)
    structureName = 'structure';
  end
  option = defopts(settings, 'Format', '');
  
  % create higher structure for structure array handling (due to function
  % subfields
  extraStruct.s = structure;
  % gain all structure subfields
  settingsSF = subfields(extraStruct);
    
  % print all subfields
  for sf = 1:length(settingsSF)
    % eval due to multiple subfields
    valueSF = eval(['extraStruct.', settingsSF{sf}]);
    switch option
      case {'value', 'values'}
        strToPrint = '';
      case {'field', 'fields'}
        strToPrint = [settingsSF{sf}(3:end), ' = '];
      otherwise
        if strcmp(structureName(end), ' ') || ~isstruct(structure)
          strToPrint = [structureName, settingsSF{sf}(3:end), ' = '];
        else
          strToPrint = [structureName, '.', settingsSF{sf}(3:end), ' = '];
        end
    end
    str = prt(str, '%s', strToPrint);
    % array settings
    if numel(valueSF) > 1 && ~ischar(valueSF)
      str = printArray(str, valueSF);
    % non-array value
    else
      str = printVal(str, valueSF);
    end
    if strcmp(option, 'value')
      if length(settingsSF) ~= 1
        str = prt(str, '\n');
      end
    else
      str = prt(str, ';\n');
    end

  end
  
  if nargout < 1
    recentlyOpened = false;
    if ischar(FID)
      % printing results to txt file
      resultname = FID;
      FID = fopen(resultname, 'w');
      assert(FID ~= -1, 'Cannot open %s !', resultname)
      recentlyOpened = true;
    end

    % print to file
      fprintf(FID, str);
      clear str

    % close opened file
    if recentlyOpened
      fclose(FID);
    end
  end
  
end

function str = prt(str, varargin)
% adds string to str
  str = [str, sprintf(varargin{:})];
end

function str = printVal(str, val)
% function checks the class of value and prints it in appropriate format

  if isempty(val)
    if iscell(val)
      str = prt(str, '{}');
    else
      str = prt(str, '[]');
    end
  elseif iscell(val) || (numel(val) > 1 && ~ischar(val))
  % cell or any kind of array (except char)
    str = printArray(str, val);
  else
    switch class(val)
      case 'char'
        str = prt(str, '''%s''', val);
      case 'double'
        % NaN and Inf verification is part of condition because they cannot 
        % be converted to logicals
        if (isnan(val) || val == Inf || (mod(val,1) && abs(val) > 1))
          str = prt(str, '%f', val);
        elseif mod(val,1)
          str = prt(str, '%g', val);
        else
          str = prt(str, '%d', val);
        end
      case 'logical'
        if val
          str = prt(str, 'true');
        else
          str = prt(str, 'false');
        end
      case 'function_handle'
        str = prt(str, '@%s', func2str(val));
      case 'struct'
        str = printStruct(str, val);
      otherwise
        str = prt(str, '%dx%d %s', size(val,1), size(val,2), class(val));
    end
  end
end

function str = printArray(str, val)
% function prints array

  % cell array
  if iscell(val)
    str = prt(str, '{');
    % first row
    str = printVal(str, val{1,1});
    for c = 2:size(val,2)
      str = prt(str, ', ');
      str = printVal(str, val{1,c});
    end
    % rest of rows
    for r = 2:size(val,1)
      str = prt(str, '; ');
      str = printVal(str, val{r,1});
      for c = 2:size(val,2)
        str = prt(str, ', ');
        str = printVal(str, val{r,c});
      end
    end
    str = prt(str, '}');
  % other arrays
  else
    str = prt(str, '[');
    % first row
    str = printVal(str, val(1,1));
    for c = 2:size(val,2)
      str = prt(str, ', ');
      str = printVal(str, val(1,c));
    end
    % rest of rows
    for r = 2:size(val,1)
      str = prt(str, '; ');
      str = printVal(str, val(r,1));
      for c = 2:size(val,2)
        str = prt(str, ', ');
        str = printVal(str, val(r,c));
      end
    end
    str = prt(str, ']');
  end
end

function str = printStruct(str, s)
% prints structure s

  sf = fieldnames(s);
  str = prt(str, 'struct(');
  str = prt(str, '''%s'', ', sf{1});
  str = printVal(str, s.(sf{1}));
  for fnum = 2 : length(sf)
    str = prt(str, ', ''%s'', ', sf{fnum});
    str = printVal(str, s.(sf{fnum}));
  end
  str = prt(str, ')');
end