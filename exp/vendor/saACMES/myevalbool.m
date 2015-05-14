% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myevalbool(s)
  if ~ischar(s) % s may not and cannot be empty
    res = s;
  else % evaluation string s
    if strncmpi(lower(s), 'yes', 3) || strncmpi(s, 'on', 2) ...
	  || strncmpi(s, 'true', 4) || strncmp(s, '1 ', 2)
      res = 1;
    elseif strncmpi(s, 'no', 2) || strncmpi(s, 'off', 3) ...
	  || strncmpi(s, 'false', 5) || strncmp(s, '0 ', 2)
      res = 0;
    else
      try res = evalin('caller', s); catch
	error(['String value "' s '" cannot be evaluated']);
      end
      try res ~= 0; catch
	error(['String value "' s '" cannot be evaluated reasonably']);
      end
    end
  end