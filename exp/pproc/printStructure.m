function printStructure(structure, FID)
% printStructure(structure, FID) prints all fields and values of 
% 'structure' to file 'FID'
%
% Input:
%   structure - structure to print
%   FID       - file identifier or filename | integer or string

  if nargin < 2
    FID = 1;
    if nargin < 1
      help printStructure
      return
    end
  end
  
  recentlyOpened = false;
  if ischar(FID)
    % printing results to txt file
    resultname = FID;
    FID = fopen(resultname, 'w');
    assert(FID ~= -1, 'Cannot open %s !', resultname)
    recentlyOpened = true;
  end
  
  % gain structure name
  structureName = inputname(1);
  if isempty(structureName)
    structureName = 'structure';
  end
  
  % create higher structure for structure array handling (due to function
  % subfields
  extraStruct.(structureName) = structure;
  % gain all structure subfields
  settingsSF = subfields(extraStruct);
    
  % print all subfields
  for sf = 1:length(settingsSF)
    % eval due to multiple subfields
    valueSF = eval(['extraStruct.', settingsSF{sf}]);
    % array settings
    if numel(valueSF) > 1 && ~ischar(valueSF)
      printArray(FID, valueSF);
    % non-array value
    else
      printVal(FID, valueSF);
    end
    if length(settingsSF) ~= 1
      fprintf(FID,'\n');
    end
  end
  
  % close opened file
  if recentlyOpened
    fclose(FID);
  end
  
end

function printVal(FID, val)
% function checks the class of value and prints it in appropriate format
% to file FID

  if isempty(val)
    if iscell(val)
      fprintf(FID,'{}');
    else
      fprintf(FID,'[]');
    end
  elseif iscell(val) || (numel(val) > 1 && ~ischar(val))
  % cell or any kind of array (except char)
    printArray(FID, val);
  else
    switch class(val)
      case 'char'
        fprintf(FID,'''%s''', val);
      case 'double'
        % NaN and Inf verification is part of condition because they cannot 
        % be converted to logicals
        if (isnan(val) || val == Inf || (mod(val,1) && abs(val) > 1))
          fprintf(FID,'%f', val);
        elseif mod(val,1)
          fprintf(FID,'%g', val);
        else
          fprintf(FID,'%d', val);
        end
      case 'logical'
        if val
          fprintf(FID,'true');
        else
          fprintf(FID,'false');
        end
      case 'function_handle'
        fprintf(FID,'@%s', func2str(val));
      case 'struct'
        printStruct(FID, val)
      otherwise
        fprintf(FID,'%dx%d %s', size(val,1), size(val,2), class(val));
    end
  end
end

function printArray(FID, val)
% function prints array to file FID

  % cell array
  if iscell(val)
    fprintf(FID,'{');
    % first row
    printVal(FID, val{1,1})
    for c = 2:size(val,2)
      fprintf(FID,', ');
      printVal(FID, val{1,c})
    end
    % rest of rows
    for r = 2:size(val,1)
      fprintf(FID,'; ');
      printVal(FID, val{r,1})
      for c = 2:size(val,2)
        fprintf(FID,', ');
        printVal(FID, val{r,c})
      end
    end
    fprintf(FID,'}');
  % other arrays
  else
    fprintf(FID,'[');
    % first row
    printVal(FID, val(1,1))
    for c = 2:size(val,2)
      fprintf(FID,', ');
      printVal(FID, val(1,c))
    end
    % rest of rows
    for r = 2:size(val,1)
      fprintf(FID,'; ');
      printVal(FID, val(r,1))
      for c = 2:size(val,2)
        fprintf(FID,', ');
        printVal(FID, val(r,c))
      end
    end
    fprintf(FID,']');
  end
end

function printStruct(FID, s)
% prints structure s to file FID

  sf = fieldnames(s);
  fprintf(FID,'struct(');
  fprintf(FID,'''%s'', ', sf{1});
  printVal(FID, s.(sf{1}))
  for fnum = 2 : length(sf)
    fprintf(FID,', ''%s'', ', sf{fnum});
    printVal(FID, s.(sf{fnum}))
  end
  fprintf(FID,')');
end

function sf = subfields(ThisStruct)
% sf = subfields(ThisStruct) returns cell array of all fields of structure 
% ThisStruct except structure names.

   sf = fieldnames(ThisStruct);
   Nsf = length(sf);
   deletesf = false(1,Nsf);
   
   for fnum = 1:Nsf
     if isstruct(ThisStruct.(sf{fnum}))
       [sRows, sCols] = size(ThisStruct.(sf{fnum}));
       % single structure case
       if sRows == 1 && sCols == 1
         sf = catSub(ThisStruct.(sf{fnum}), sf, fnum, '');
       % structure vector
       elseif sRows == 1 || sCols == 1
         longerDimVal = max(sRows, sCols);
         for v = 1:longerDimVal
           sf = catSub(ThisStruct.(sf{fnum})(v), sf, fnum, ['(', num2str(v), ')']);
         end
       % structure array
       else
         for r = 1:sRows
           for c = 1: sCols
             sf = catSub(ThisStruct.(sf{fnum})(r,c), sf, fnum, ['(', num2str(r), ',', num2str(c), ')']);
           end
         end
       end
       % mark field as possible to delete
       deletesf(fnum) = true;
     end
   end
   % delete higher structure names
   sf(deletesf) = [];

end

function sf = catSub(ThisStruct, sf, fnum, structStr)
% function concatenates structure subfields
  cn = subfields(ThisStruct);
  sf = cat(1, sf, strcat(sf{fnum}, structStr, '.', cn));
end