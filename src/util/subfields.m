function sf = subfields(thisStruct)
% sf = subfields(ThisStruct) returns cell array of all field and subfield
% names of structure ThisStruct except structure names.
%
% See Also:
%   printStructure

   sf = fieldnames(thisStruct);
   Nsf = length(sf);
   deletesf = false(1,Nsf);
   
   for fnum = 1:Nsf
     if isstruct(thisStruct.(sf{fnum}))
       [sRows, sCols] = size(thisStruct.(sf{fnum}));
       % single structure case
       if sRows == 1 && sCols == 1
         sf = catSub(thisStruct.(sf{fnum}), sf, fnum, '');
       % structure vector
       elseif sRows == 1 || sCols == 1
         longerDimVal = max(sRows, sCols);
         for v = 1:longerDimVal
           sf = catSub(thisStruct.(sf{fnum})(v), sf, fnum, ['(', num2str(v), ')']);
         end
       % structure array
       else
         for r = 1:sRows
           for c = 1: sCols
             sf = catSub(thisStruct.(sf{fnum})(r,c), sf, fnum, ['(', num2str(r), ',', num2str(c), ')']);
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