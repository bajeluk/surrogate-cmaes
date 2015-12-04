mailTo = 'bajeluk@gmail.com,z.pitra@gmail.com';
sendMail = false;

% Divide instances to machines
params = [bbParamDef, sgParamDef, cmParamDef];
nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
nMachines = length(machines);

% Estimate the running time based on the dimensions and inline function
% provided to divideTasksForMachines()
dimensions = zeros(1,nCombinations);
models = zeros(1,nCombinations);
estTimes = zeros(1,nCombinations);
for id = 1:nCombinations
  [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
  dimensions(id) = bbParams.dimensions;
  models(id) = isfield(sgParams, 'modelType') && strcmpi(sgParams.modelType, 'rf');
end
estTimes = dimensions; % + 3*models.*dimensions;
cellCombsForMachines = divideTasksForMachines(nMachines, estTimes, @(x) x.^(1.3));

% Generate .sh scripts
fNameMng = [exp_id '_manager.sh'];
fMng = fopen([exppath filesep fNameMng], 'w'); 
fprintf(fMng, '#!/bin/sh\n');
fprintf(fMng, '# Manager for experiment "%s", created on %s\n', exp_id, datestr(now,'yyyy-mm-dd HH:MM:SS'));

for i = 1:min([nMachines nCombinations])
  machine = machines{i};
  fName = [exppath filesep exp_id '_' machine '.sh'];

  % Logging of Matlab output
  if (exist('logMatlabOutput', 'var') && logMatlabOutput)
    logFile = ['$EXPPATH' filesep exp_id '__log__' machine '.txt'];
    logString = ' | tee -a $LOGS';
  else
    logFile = '/dev/null';
    logString = '';
  end

  fprintf(fMng, 'ssh %s@%s "screen -d -m \\"%s\\""\n', login, machine, fName);

  fid = fopen(fName, 'w'); 
  fprintf(fid, '#!/bin/sh\n');
  fprintf(fid, '# Script for experiment "%s" on machine "%s", created on %s\n', exp_id, machine, datestr(now,'yyyy-mm-dd HH:MM:SS'));
  fprintf(fid, 'cd ~/prg/surrogate-cmaes; ulimit -t unlimited;\n');
  fprintf(fid, 'EXPPATH_SHORT="%s"\n', exppath_short);
  fprintf(fid, 'EXPPATH="$EXPPATH_SHORT%s"\n', [filesep exp_id]);
  fprintf(fid, 'LOGS="%s"\n\n', logFile);
  fprintf(fid, 'testMatlabFinished () {\n');
  fprintf(fid, '  if [ -f "$1" -a "$1" -nt "$EXPPATH%s" ]; then\n', [filesep fNameMng]);
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` " " **%s** at [%s] $2 / $3 succeeded. >> ~/WWW/phd/cmaes.txt\n', exp_id, machine);
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` " " **%s** at [%s] $2 / $3 succeeded. >> %s\n', exp_id, machine, logFile);
  fprintf(fid, '  else\n');
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` " " **%s** at [%s] $2 / $3 !!!  FAILED !!! >> ~/WWW/phd/cmaes.txt\n', exp_id, machine);
  fprintf(fid, '    echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` " " **%s** at [%s] $2 / $3 !!!  FAILED !!! >> %s\n', exp_id, machine, logFile);
  % Sending mail on error
  if (sendMail)
    fprintf(fid, '    mail -s "[surrogate-cmaes] Script failed =%s= [%s] $2 / $3" "%s" <<EOM\n', exp_id, machine, mailTo);
    fprintf(fid, 'Script **%s** has failed on machine [%s].\n', exp_id, machine);
    fprintf(fid, '`date "+%%Y-%%m-%%d %%H:%%M:%%S"`\n');
    fprintf(fid, 'EOM\n');
  end
  fprintf(fid, '  fi\n');
  fprintf(fid, '}\n');
  idFrom1 = 0;
  for id = cellCombsForMachines{i}
    idFrom1 = idFrom1 + 1;
    nCombsForThisMachine = length(cellCombsForMachines{i});
    [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);

    lastResultsFileID = [num2str(bbParams.functions(end)) '_' num2str(bbParams.dimensions(end)) 'D_' num2str(id)];
    resultsFile = ['$EXPPATH' filesep exp_id '_results_' lastResultsFileID '.mat'];

    fprintf(fid, '\necho "###########################################"%s\n', logString);
    fprintf(fid, 'echo "     Matlab call %d / %d"%s\n', idFrom1, nCombsForThisMachine, logString);
    fprintf(fid, 'echo ""%s\n', logString);
    fprintf(fid, 'echo "  dim(s): %s    f(s): %s    N(inst): %d"%s\n', num2str(bbParams.dimensions), num2str(bbParams.functions), length(bbParams.instances), logString);
    if (isfield(sgParams, 'modelType'))
      model = sgParams.modelType;
    else
      model = 'NONE';
    end
    fprintf(fid, 'echo "  model: %s    maxfunevals: %s"%s\n', model, bbParams.maxfunevals, logString);
    fprintf(fid, 'echo "###########################################"%s\n', logString);
    fprintf(fid, 'nice -n 19 %s -nodisplay -r "bbob_test_01(%d, ''%s'', ''$EXPPATH_SHORT''); exit(0);" 2>&1 %s\n', matlabcommand, id, exp_id, logString);
    fprintf(fid, 'testMatlabFinished "%s" %d %d\n', resultsFile, idFrom1, nCombsForThisMachine);
  end
  
  fprintf(fid, 'echo `date "+%%Y-%%m-%%d %%H:%%M:%%S"` " " **%s** at [%s] ==== FINISHED ==== >> ~/WWW/phd/cmaes.txt\n', exp_id, machine);

  fclose(fid);
  if isunix
    fileattrib(fName, '+x');
  end
end

fclose(fMng);
if isunix
  fileattrib([exppath filesep fNameMng], '+x');
end
