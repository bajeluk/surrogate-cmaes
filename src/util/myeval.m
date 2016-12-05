% MYEVAL - evaluate or use the expression
% author: N. Hansen
%
function res = myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
