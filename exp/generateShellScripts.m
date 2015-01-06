% Divide instances to machines
params = [bbParamDef, sgParamDef, cmParamDef];
nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
nMachines = length(machines);

combsPerMachine = ceil(nCombinations / nMachines);
startIdxs = 1:combsPerMachine:(nCombinations-1);
endIdxs =   combsPerMachine:combsPerMachine:nCombinations;

% Generate .sh scripts
for i = 1:nMachines
  machine = machines{i};
  fid = fopen([exppath filesep exp_id '_' machine '.sh'], 'w'); 
  fprintf(fid, '#!/bin/sh\n');
  for id = startIdxs(i):endIdxs(i)
    matlabCall = ['ssh ' login '@' machine ' "cd ~/prg/surrogate-cmaes; ulimit -t unlimited; nice -n 19 ' matlabcommand ' -r \\"bbob_test_01(' num2str(id) ', ''' exp_id ''', ''' exppath_short '''); exit\\""'];
    fprintf(fid, [matlabCall '\n']);
    disp(matlabCall);
  end
  fclose(fid);
end

