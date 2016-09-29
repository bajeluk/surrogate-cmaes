%READSCMAESTXTLOG Reads S-CMA-ES log from .txt into Matlab table
function tab = readScmaesTxtLog(exp_id, fun, dim, id)
  cwd = fileparts(mfilename('fullpath'));
  path = [cwd '/../experiments/' exp_id '/bbob_output/'];
  filename = sprintf('%s/%s_log_%d_%dD_%d.dat', path, exp_id, fun, dim, id);
  randNum = 1000+randi(1000);
  tmpFile = ['/tmp/scmaes_dat_' num2str(randNum)];
  cmd = ['sed -n ''2,/^#  f/p'' ' filename ' | sed ''s/^ *//;s/  */ /g;1s/^# //;s/ *$//;$d'' > ' tmpFile];
  system(cmd);
  tab = readtable(tmpFile, 'Delimiter', ' ', 'ReadVariableNames', true);
  delete(tmpFile);
end
