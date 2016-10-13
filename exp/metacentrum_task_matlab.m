function metacentrum_task_matlab(exp_id, exppath_short, id, varargin)
% metacentrum_task_matlab(exp_id, exppath_short, id, varargin) calls the
% individual task (with the chosen task_id). 
%
% Important: Fuction expInit(exp_id) has to be called before the 
%            metacentrum_task_matlab function to initialize the experiment.
%
% Input:
%   exp_id        - experiment ID (exp_id) (string)
%   exppath_short - absolute path to the directory exp/experiments (string)
%   id            - task_id of the task to be run (integer)
%                 - specifies parameter-values being used for particular 
%                   run
%   varargin      - optional argument specifying additional settings 
%                   (structure array)
%
% See Also:
%   metacentrum_master_template

  if nargin < 3
    help metacentrum_task_matlab
    return
  end

  % NFS file for logging results
  USE_FILELOG = 0;
  LOGFILENAME='log.txt';
  DEFAULTLOGFILE = '/storage/plzen1/home/bajeluk/public/log.txt';

  % FTP connection for logging results, NFS file log is needed for FTP!
  USE_FTPLOG = 1  &&  USE_FILELOG;
  FTP_HOSTNAME='optim.wz.cz';
  FTP_USERNAME='optim.wz.cz';
  FTP_PASS='metacentrum';

  % *__log__* file -- set empty to suppress these outputs
  FILESTDOUT = [];
  % FILESTDOUT = [exppath_short filesep exp_id filesep exp_id '__log__' num2str(id) '.txt'];

  % all params are strings if compiled and called from shell
  if (ischar(id)) id = str2num(id); end

  % rename queued_file to computing_file
  queuedFile = [exppath_short filesep exp_id filesep 'queued_' num2str(id)];
  calculatingFile = [exppath_short filesep exp_id filesep 'calculating_' num2str(id)];
  if (exist(queuedFile, 'file'))
    movefile(queuedFile, calculatingFile);
  end

  % parameters of the current job
  load([exppath_short filesep exp_id filesep 'scmaes_params.mat'], 'bbParamDef', 'sgParamDef', 'cmParamDef');
  [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  fun = bbParams.functions(end);
  dim = bbParams.dimensions(end);
  if (isfield(sgParams, 'modelType'))
    model = sgParams.modelType;
  else
    model = 'NONE';
  end

  % setting up paths
  EXPPATH = [exppath_short filesep exp_id];
  RESULTSFILE = [EXPPATH '/' exp_id '_results_' num2str(fun) '_' num2str(dim) 'D_' num2str(id) '.mat'];
  if (isempty(getenv('SCRATCHDIR')))
    OUTPUTDIR = [];     % set OUTPUTDIR empty if $SCRATCHDIR var does not exist
  else
    OUTPUTDIR = [getenv('SCRATCHDIR') filesep 'job_output'];
  end

  % metaOpts -- structure with info about current Task and Metacentrum environ. variables
  metaOpts.logdir = '';
  metaOpts.machine = '';
  metaOpts.model = model;
  if (nargin >= 4)
    for fname = fieldnames(varargin{1})'
      metaOpts.(fname{1}) = varargin{1}.(fname{1});
    end
  end
  nodeFile = fopen(getenv('PBS_NODEFILE'), 'r');
  if (nodeFile > 0)
    metaOpts.machine = fgetl(nodeFile);
    fclose(nodeFile);
  end
  metaOpts.nInstances = length(bbParams.instances);

  if (USE_FILELOG)
    if (~ isfield(metaOpts, 'logdir') || isempty(metaOpts.logdir))
      LOGFILE = DEFAULTLOGFILE;
      metaOpts.logdir = '/storage/plzen1/home/bajeluk/public/';
    else
      LOGFILE = [metaOpts.logdir filesep LOGFILENAME];
    end
  end

  % # relict from shell script
  % # clean up the lock-file
  % trap "rm -f $EXPPATH/queued_$ID" TERM EXIT

  % This CD is now wrong, Matlab should be called from a SCRACHDIR
  % cd([exppath_short filesep '..' filesep '..']);
  startup;

  % STDOUT logging
  if (~isempty(FILESTDOUT))
    fout = fopen(FILESTDOUT, 'a');
    fprintf(fout, '###########################################\n');
    fprintf(fout, '     Matlab call id=%d\n\n', id);
    fprintf(fout, '  dim(s): %d    f(s): %d    N(inst): %d\n', dim, fun, metaOpts.nInstances);
    fprintf(fout, '  model: %s\n', metaOpts.model);
    fprintf(fout, '###########################################\n');
    fclose(fout);
    % TODO: forward all the text output into a file, but it has to be on a local SCRATCH
  end

  datest = datestr(now,'yyyy-mm-dd HH:MM:ss');
  % LOGFILE loggig start and success of each experiment ID into local/NFS file
  if (USE_FILELOG)
    flog = fopen(LOGFILE, 'a');
    fprintf(flog, '%s  **%s** at [%s] %d started.\n', datest, exp_id, metaOpts.machine, id);
    fclose(flog);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CORE COMPUTATION (begin)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bbob_test_01(id, exp_id, exppath_short, OUTPUTDIR);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CORE COMPUTATION (end)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  % STDOUT logging
  if (~isempty(FILESTDOUT))
    fout = fopen(FILESTDOUT, 'a');
    if (exist(RESULTSFILE, 'file'))
      fprintf(fout, '%s  **%s** at [%s] %d succeeded.\n', datest, exp_id, metaOpts.machine, id);
    else
      fprintf(fout, '%s  **%s** at [%s] %d !!!  FAILED  !!!\n', datest, exp_id, metaOpts.machine, id);
    end
    fclose(fout);
  end

  % LOGFILE loggig start and success of each experiment ID into local/NFS file
  if (USE_FILELOG)
    flog = fopen(LOGFILE, 'a');
    if (exist(RESULTSFILE, 'file'))
      % if [ "$3" -eq 0  -a  -f "$1"  -a  "$1" -nt "$FILEMANAGER" ]; then
      fprintf(flog, '%s  **%s** at [%s] %d succeeded.\n', datest, exp_id, metaOpts.machine, id);
    else
      fprintf(flog, '%s  **%s** at [%s] %d !!!  FAILED  !!!\n', datest, exp_id, metaOpts.machine, id);
      % # tail -n 60 "$FILESTDOUT" | mail -s "Metacentrum: chyba v uloze $PBS_JOBID $PBS_JOBNAME" $USER@arien.metacentrum.cz
    end
    fprintf(flog, '%s  **%s** at [%s] ==== FINISHED ====\n', datest, exp_id, metaOpts.machine);
    fclose(flog);
  end

  % copy to LOGFILE via FTP
  if (USE_FTPLOG)
    ftpstr = sprintf('ftp -n %s <<EOD\nuser %s %s\nlcd %s\nput %s\nbye\nEOD\n', FTP_HOSTNAME, FTP_USERNAME, FTP_PASS, metaOpts.logdir, LOGFILENAME);
    system(ftpstr);
  end

  % copy the BBOB results onto persistant storage if outside EXPPATH
  % ( should be done already in bbob_test_01() )
  if (~isempty(OUTPUTDIR) && ~strcmpi(OUTPUTDIR, EXPPATH) && isunix)
    % copy the output to the final storage (if OUTPUTDIR and EXPPATH differs)
    system(['cp -pR ' OUTPUTDIR '/* ' EXPPATH '/']);
  end

end
