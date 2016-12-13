%READSCMAESTXTLOG Reads S-CMA-ES log from .txt into Matlab table
%
% returns data in the Matlab table
%
% optional argument: the number of instance to load from the text log file
%   (its order, not the exact number of the instance)
function tab = readScmaesTxtLog(exp_id, fun, dim, id, varargin)
  cwd = fileparts(mfilename('fullpath'));
  path = [cwd '/../experiments/' exp_id '/bbob_output'];
  filename = sprintf('%s/%s_log_%d_%dD_%d.dat', path, exp_id, fun, dim, id);
  randNum = 1000+randi(1000);
  tmpFile = ['/tmp/scmaes_dat_' num2str(randNum)];

  sed_deleteAll = '';
  if (~isempty(varargin) && ~isempty(varargin{1}) && isnumeric(varargin{1}) && varargin{1} > 0 && varargin{1} < 16)
    instance = varargin{1};
    sed_delete_instances_string = '1,/^#  f/d;';
    for i = 1:(instance-1)
      sed_deleteAll = [sed_deleteAll sed_delete_instances_string];
    end
  end

  cmd = ['sed -n ''1d;' sed_deleteAll '1,/^#  f/p'' ' filename ' | sed ''s/^ *//;s/  */ /g;1s/^# //;s/ *$//;$d'' > ' tmpFile];
  system(cmd);
  tab = readtable(tmpFile, 'Delimiter', ' ', 'ReadVariableNames', true);
  delete(tmpFile);
end
