% TODO:
%  - oddelit SRC/EXPERIMENT-DEFINITION a OUTPUT directory

fNameMetacentrumTaskTemplate = ['exp' filesep 'metacentrum_task_template.sh'];
fNameMetacentrumTaskTemplateBinary = ['exp' filesep 'metacentrum_binary_task_template.sh'];
fNameTask = [exppath filesep exp_id '_task.sh'];
fNameTaskBinary = [exppath filesep 'binary_task.sh'];
fNameTaskShell = ['$EXPPATH_SHORT/$EXPID/${EXPID}_task.sh'];
fNameMng =  [exppath filesep exp_id '_manager.sh'];
fNameAllIds = [exppath filesep 'allids.txt'];
logDir = '/storage/plzen1/home/bajeluk/public';

% walltime = '1d';
% nCPU = 1;
% ram = '2gb';
% scratch = '1gb';
% matlabToolboxes = ',matlab_Statistics_Toolbox=1,matlab_Optimization_Toolbox=1';
% matlabParams = '-singleCompThread -nodisplay -nodesktop';
% matlabCall = 'matlab';
% % logFile = '$EXPPATH_SHORT/../../log/experiments.txt';

% Divide instances to machines
params = [bbParamDef, sgParamDef, cmParamDef];
nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);

% Export textfile with all the combinations' IDs
fIds = fopen(fNameAllIds, 'w');
if (fIds > 0)
  for id = 1:nCombinations
    fprintf(fIds, '%d ', id);
  end
  fprintf(fIds, '\n');
end

% Estimate the running time based on the dimensions and inline function
% provided to divideTasksForMachines()
dimensions = zeros(1,nCombinations);
models = zeros(1,nCombinations);
for id = 1:nCombinations
  [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  dimensions(id) = bbParams.dimensions;
  models(id) = isfield(sgParams, 'modelType') && strcmpi(sgParams.modelType, 'rf');
end
estTimes = dimensions; % + 3*models.*dimensions;

% nMachines = length(machines);
% cellCombsForMachines = divideTasksForMachines(nMachines, estTimes, @(x) x.^(1.3));

% copy task shell script into the directory with the experiment
copyfile([pathstr filesep '..' filesep '..' filesep fNameMetacentrumTaskTemplate], fNameTask);
copyfile([pathstr filesep '..' filesep '..' filesep fNameMetacentrumTaskTemplateBinary], fNameTaskBinary);
if (isunix)
  fileattrib(fNameTask, '+x');
  fileattrib(fNameTaskBinary, '+x');
end

% Generate managing .sh script
fMng = fopen(fNameMng, 'w'); 
fprintf(fMng, '#!/bin/sh\n');
fprintf(fMng, '# Manager for experiment "%s", created on %s\n', exp_id, datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fMng, '\n');
fprintf(fMng, 'export EXPID\n');
fprintf(fMng, 'export ID\n');
fprintf(fMng, 'export FUN\n');
fprintf(fMng, 'export DIM\n');
fprintf(fMng, 'export NINSTANCES\n');
fprintf(fMng, 'export MODEL\n');
fprintf(fMng, 'export EXPPATH_SHORT\n');
fprintf(fMng, 'export LOGDIR\n');
fprintf(fMng, '\n');
fprintf(fMng, 'EXPID="%s"\n', exp_id);
fprintf(fMng, 'EXPPATH_SHORT="%s"\n', exppath_short);
fprintf(fMng, 'LOGDIR="%s"\n', logDir);
fprintf(fMng, '\n');

for id = 1:nCombinations
  [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  fileID = [num2str(bbParams.functions(end)) '_' num2str(bbParams.dimensions(end)) 'D_' num2str(id)];

  if (isfield(sgParams, 'modelType'))
    model = sgParams.modelType;
  else
    model = 'NONE';
  end

  % set the parameters for the task
  fprintf(fMng, 'ID=%d\n', id);
  fprintf(fMng, 'FUN=%d; DIM=%d; NINSTANCES=%d; MODEL="%s"\n', bbParams.functions(end), bbParams.dimensions(end), length(bbParams.instances), model);
  % add queueing to the "manager"
  fprintf(fMng, 'qsub -N "${EXPID}_${ID}" -v EXPID,ID,FUN,DIM,NINSTANCES,MODEL,EXPPATH_SHORT,LOGDIR %s\n', fNameTaskShell);
  fprintf(fMng, 'if [ ! $? -eq 0 ] ; then\n');
  fprintf(fMng, '  echo "Nepodarilo se zadat ulohu segment ${ID}! Koncim."; exit 1\n');
  fprintf(fMng, 'else\n');
  fprintf(fMng, '  touch "$EXPPATH_SHORT/$EXPID/queued_$ID"\n');
  fprintf(fMng, 'fi\n');
  fprintf(fMng, '\n');

  %{
  % create the file for the script itself
  fid = fopen(fNameTask, 'w'); 
  fprintf(fid, '#!/bin/sh\n');
  fprintf(fid, '#PBS -N %s_%s\n', exp_id, fileID);
  fprintf(fid, '#PBS -l walltime=%s', walltime);
  fprintf(fid, '#PBS -l nodes=1:ppn=%d\n', nCPU);
  fprintf(fid, '#PBS -l mem=%s\n', ram);
  fprintf(fid, '#PBS -l scratch=%s\n', scratch);
  fprintf(fid, '#PBS -l matlab=1%s\n', matlabToolboxes);
  % join the STDID and STDERR
  % fprintf(fid, '#PBS -j oe\n'
  % and send it by e-mail
  % fprintf(fid, '#PBS -m e\n'

  fprintf(fid, '#\n');
  fprintf(fid, '# Script for experiment "%s", created on %s\n', exp_id, datestr(now,'yyyy-mm-dd HH:MM:SS'));
  fprintf(fid, '#\n');
  fprintf(fid, 'EXPID="%s"\n', exp_id);
  fprintf(fid, 'ID=%d\n', id);
  fprintf(fid, 'FUN=%d\n', bbParams.functions(end));
  fprintf(fid, 'DIM=%d\n', bbParams.dimensions(end));
  fprintf(fid, 'NINSTANCES=%d\n', length(bbParams.instances));
  fprintf(fid, 'MODEL="%s"\n', model);
  fprintf(fid, 'EXPPATH_SHORT="%s"\n', exppath_short);
  fprintf(fid, 'LOGFILE="%s"\n', logFile);
  fprintf(fid, '\n');
  fprintf(fid, 'HOSTNAME=optim.wz.cz\n');
  fprintf(fid, 'USERNAME=optim.wz.cz\n');
  fprintf(fid, 'PASS=metacentrum\n');
  fprintf(fid, 'EXPPATH="$EXPPATH_SHORT/$EXPID"\n');
  fprintf(fid, 'OUTPUTDIR="$EXPPATH"\n');
  fprintf(fid, 'RESULTSFILE="${EXPID}_results_${FUN}_${DIM}D_${ID}.mat"\n');
  fprintf(fid, 'FILESTDOUT="$OUTPUTDIR/${EXPID}__log__${ID}.txt"\n');
  fprintf(fid, 'FILEMANAGER="${EXPPATH}/${EXPID}_manager.sh"\n');
  fprintf(fid, 'MATLABCALL="%s"\n', matlabCall);
  fprintf(fid, 'MATLABPARAMS="%s"\n', matlabParams);
  fprintf(fid, 'MACHINE=`head -1 $PBS_NODEFILE`\n\n');
  fprintf(fid, 'testMatlabFinished () {\n');
  fprintf(fid, '  if [ "$2" -eq 0  -a  -f "$1"  -a  "$1" -nt "$FILEMANAGER" ]; then\n');
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` "  **$EXPID** at [$MACHINE] $2 succeeded." >> $LOGFILE\n');
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` "  **$EXPID** at [$MACHINE] $2 succeeded." >> $FILESTDOUT\n');
  fprintf(fid, '  else\n');
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` "  **$EXPID** at [$MACHINE] $2 !!!  FAILED !!!" >> $LOGFILE\n');
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` "  **$EXPID** at [$MACHINE] $2 !!!  FAILED !!!" >> $FILESTDOUT\n');
  fprintf(fid, '  fi\n');
  fprintf(fid, '  ftp -n $HOSTNAME <<EOD\n');
  fprintf(fid, 'user $USERNAME $PASS\n');
  fprintf(fid, 'put $LOGFILE\n');
  fprintf(fid, 'bye\n');
  fprintf(fid, 'EOD\n');
  fprintf(fid, '}\n\n');

  fprintf(fid, 'cd "$EXPPATH_SHORT/../.."; ulimit -t unlimited;\n');

  fprintf(fid, '\necho "###########################################"%s\n', logString);
  fprintf(fid, 'echo "     Matlab call id=${ID}"%s\n', logString);
  fprintf(fid, 'echo ""%s\n', logString);
  fprintf(fid, 'echo "  dim(s): $DIM    f(s): $FUN    N(inst): $NINSTANCES"%s\n', logString);
  fprintf(fid, 'echo "  model: $MODEL"%s\n', logString);
  fprintf(fid, 'echo "###########################################"%s\n', logString);
  fprintf(fid, 'echo "exit(1)" | $MATLABCALL $MATLABPARAMS -r "bbob_test_01($ID, ''$EXPID'', ''$EXPPATH_SHORT'', ''$OUTPUTDIR''), exit(0)" %s 2>&1\n', logString);
  fprintf(fid, '\n');
  fprintf(fid, 'testMatlabFinished "$OUTPUTDIR/$RESULTSFILE" $ID $?\n');
  
  fprintf(fid, '\necho `date "+%%Y-%%m-%%d %%H:%%M:%%S"` "  **$EXPID** at [$MACHINE] ==== FINISHED ====" >> $LOGFILE\n');

  fclose(fid);
  if (isunix) fileattrib(fNameTask, '+x'); end
  %}
end

fclose(fMng);
if (isunix) fileattrib(fNameMng, '+x'); end
