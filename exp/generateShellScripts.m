% Divide instances to machines
params = [bbParamDef, sgParamDef, cmParamDef];
nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
nMachines = length(machines);

combsPerMachine = ceil(nCombinations / nMachines);
startIdxs = 1:combsPerMachine:(nCombinations-1);
endIdxs =   combsPerMachine:combsPerMachine:nCombinations;

% Generate .sh scripts
fNameMng = [exppath filesep exp_id '_manager.sh'];
fMng = fopen(fNameMng, 'w'); 
fprintf(fMng, '#!/bin/sh\n');
fprintf(fMng, '# Manager for experiment "%s", created on %s\n', exp_id, char(datetime('now')));

for i = 1:nMachines
  machine = machines{i};
  fName = [exppath filesep exp_id '_' machine '.sh'];

  fprintf(fMng, 'ssh %s@%s "screen -d -m \\"%s\\""\n', login, machine, fName);

  fid = fopen(fName, 'w'); 
  fprintf(fid, '#!/bin/sh\n');
  fprintf(fid, '# Script for experiment "%s" on machine "%s", created on %s\n', exp_id, machine, char(datetime('now')));
  fprintf(fid, 'cd ~/prg/surrogate-cmaes; ulimit -t unlimited;\n');
  for id = startIdxs(i):endIdxs(i)
    fprintf(fid, 'echo ###########################################\n');
    fprintf(fid, 'echo      Matlab call %d / %d\n', id - startIdxs(i)+1, combsPerMachine);
    fprintf(fid, 'echo ###########################################\n');
    matlabCall = ['nice -n 19 ' matlabcommand ' -nodisplay -r "bbob_test_01(' num2str(id) ', ''' exp_id ''', ''' exppath_short '''); exit"'];
    fprintf(fid, [matlabCall '\n']);
    disp(matlabCall);
  end
  fclose(fid);
  fileattrib(fName, '+x');
end

fclose(fMng);
fileattrib(fNameMng, '+x');
