function metacentrum_task_matlab(exp_id, exppath_short, fun, dim, id, metaOpts)

  % metaOpts.logdir = '';
  % metaOpts.model = '';
  % metaOpts.nInstances = '';
  [~, metaOpts.machine] = system('head -1 $PBS_NODEFILE');
  metaOpts.machine = metaOpts.machine(1:end-1);

  FTP_HOSTNAME='optim.wz.cz';
  FTP_USERNAME='optim.wz.cz';
  FTP_PASS='metacentrum';
  LOGFILENAME='log.txt';
  DEFAULTLOGFILE = '/storage/plzen1/home/tosovvoj/public/log.txt';
  EXPPATH = [exppath_short filesep exp_id];
  OUTPUTDIR = getenv('SCRATCHDIR');
  RESULTSFILE = [EXPPATH '/' exp_id '_results_' num2str(fun) '_' num2str(dim) 'D_' num2str(id) '.mat'];
  % this could be changed to [OUTPUTDIR '/' ...]:
  FILESTDOUT = [EXPPATH '/' exp_id '__log__' num2str(id) '.txt'];
  % FILEMANAGER="${EXPPATH}/${EXPID}_manager.sh"
  % MATLABCALL="matlab"
  % MATLABPARAMS="-singleCompThread -nodisplay -nodesktop"
  % MACHINE=`head -1 $PBS_NODEFILE`

  if (~ isfield(metaOpts, 'logdir') || isempty(metaOpts.logdir))
    LOGFILE = DEFAULTLOGFILE;
    metaOpts.logdir = '/storage/plzen1/home/tosovvoj/public/';
  else
    LOGFILE = [metaOpts.logdir filesep LOGFILENAME];
  end

  % # clean up the lock-file
  % trap "rm -f $EXPPATH/queued_$ID" TERM EXIT

  cd([exppath_short filesep '..' filesep '..']);
  startup;

  fout = fopen(FILESTDOUT, 'a');

  fprintf(fout, '###########################################\n');
  fprintf(fout, '     Matlab call id=%d\n\n', id);
  fprintf(fout, '  dim(s): %d    f(s): %d    N(inst): %d\n', dim, fun, metaOpts.nInstances);
  fprintf(fout, '  model: %s\n', metaOpts.model);
  fprintf(fout, '###########################################\n');

  datest = datestr(now,'yyyy-mm-dd HH:MM:ss');
  flog = fopen(LOGFILE, 'a');
  fprintf(flog, '%s  **%s** at [%s] %d started.\n', datest, exp_id, metaOpts.machine, id);
  fclose(fout);
  fclose(flog);

  bbob_test_01(id, exp_id, exppath_short, OUTPUTDIR);

  fout = fopen(FILESTDOUT, 'a');
  flog = fopen(LOGFILE, 'a');
  datest = datestr(now,'yyyy-mm-dd HH:MM:ss');
  if (exist(RESULTSFILE, 'file'))
    % if [ "$3" -eq 0  -a  -f "$1"  -a  "$1" -nt "$FILEMANAGER" ]; then
    fprintf(flog, '%s  **%s** at [%s] %d succeeded.\n', datest, exp_id, metaOpts.machine, id);
    fprintf(fout, '%s  **%s** at [%s] %d succeeded.\n', datest, exp_id, metaOpts.machine, id);
  else
    fprintf(flog, '%s  **%s** at [%s] %d !!!  FAILED  !!!\n', datest, exp_id, metaOpts.machine, id);
    fprintf(fout, '%s  **%s** at [%s] %d !!!  FAILED  !!!\n', datest, exp_id, metaOpts.machine, id);
    % # tail -n 60 "$FILESTDOUT" | mail -s "Metacentrum: chyba v uloze $PBS_JOBID $PBS_JOBNAME" $USER@arien.metacentrum.cz
  end
  fprintf(flog, '%s  **%s** at [%s] ==== FINISHED ====\n', datest, exp_id, metaOpts.machine);
  fclose(flog);
  fclose(fout);

  ftpstr = sprintf('ftp -n %s <<EOD\nuser %s %s\nlcd %s\nput %s\nbye\nEOD\n', FTP_HOSTNAME, FTP_USERNAME, FTP_PASS, metaOpts.logdir, LOGFILENAME);
  system(ftpstr);
  if (~strcmpi(OUTPUTDIR, EXPPATH))
    % copy the output to the final storage (if OUTPUTDIR and EXPPATH differs)
    system(['cp -r ' OUTPUTDIR '/* ' EXPPATH '/']);
  end
end
