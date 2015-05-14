% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = isoctave
% any hack to find out whether we are running octave
  s = version;
  res = 0;
  if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
  end