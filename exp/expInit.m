function expInit(exp_id)
% EXPINIT -- S-CMA-ES experiment initialization
%
% EXP_ID -- string with an unique name of the experiment
%
% expInit('exp_name') has to be run before running the experiment.
% The experiment definitions are located in m-files 'exp/experiments/[exp_id].m'.
% New syntax of these experiment definitions is rather straightforward and
% this function EXPINIT prepares the experiment-definition-data in the same
% format as the old version of experiment definitions and script
% 'generateShellScriptsMetacentrum'
%
% SEE ALSO: generateStructOpts.m, getParamIndexVector.m, expUtilTest.m

  exppath_short = fileparts(mfilename('fullpath'));

  expFile = fullfile(exppath_short, 'experiments', [exp_id '.m']);
  if (~exist(expFile, 'file'))
    error('Error: experiment file "%s" does not exist.', expFile);
  end

  % load the core experiment configurations
  run(expFile);

  % BBOB/COCO settings
  bbParamDef = generateParamDefStruct(bbobParams);

  % Surrogate-modelling settings
  sgParamDef = generateParamDefStruct(surrogateParams);
  % prepare a directory for experiments settings and results: exp/experiments/[exp_id]
  experimentPath = fullfile(exppath_short, 'experiments', exp_id);
  mkdir(experimentPath);
  % save this directory as a new field to sgParamDef
  sgParamDef(end+1).name = 'experimentPath';
  sgParamDef(end).values = { experimentPath };

  % Model settings (is historically part of sgParamDef)
  % -- enumerate all the possible combinations
  sgParamDef(end+1).name = 'modelOpts';
  multiStruct = generateStructOpts(modelParams);
  sgParamDef(end).values = convertMultiStructToCell(multiStruct);

  % CMA-ES settings
  cmParamDef = generateParamDefStruct(cmaesParams);

  % save the definition of the experiment into 'scmaes_params.mat'
  save([experimentPath filesep 'scmaes_params.mat'], ...
    'bbParamDef', 'sgParamDef', 'cmParamDef', 'exp_id', 'exppath_short', 'exp_description', 'logDir');

  % prepare other directories and files in the 'experimentPath'
  [s,mess,messid] = mkdir([experimentPath filesep 'cmaes_results']);
  params = [bbParamDef, sgParamDef, cmParamDef];
  nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
  prepareMetacentrumFiles(exppath_short, experimentPath, exp_id, nCombinations);
end


function defs = generateParamDefStruct(carray)
  assert(mod(length(carray), 2) == 0, 'The number of elements in field-value cell array is not even!')
  defs = struct();
  for i = 1:(length(carray)/2)
    defs(i).name   = carray{(i-1)*2 + 1};
    defs(i).values = carray{i*2};
  end
end


function c = convertMultiStructToCell(st)
  c = cell(1,length(st));
  for i = 1:length(st)
    c{i} = st(i);
  end
end


function prepareMetacentrumFiles(expDir, experimentPath, exp_id, nCombinations)
  % expDir         = path to   'exp/'
  % experimentPath = path to   'exp/experiment/exp_id'
  % exp_id         = string with the unique experiment identifier

  % copy MCR-compiled task shell script into the directory with the experiment
  fNameMetacentrumTaskTemplateBinary = [expDir filesep 'metacentrum_binary_task_template.sh'];
  fNameTaskBinary   = [experimentPath filesep 'binary_task.sh'];
  copyfile(fNameMetacentrumTaskTemplateBinary,  fNameTaskBinary);
  if (isunix)
    fileattrib(fNameTaskBinary, '+x');
  end

  % Export a textfile 'allids.txt' with all the combinations' IDs
  fNameAllIds       = [experimentPath filesep 'allids.txt'];
  fIds = fopen(fNameAllIds, 'w');
  if (fIds > 0)
    for id = 1:nCombinations
      fprintf(fIds, '%d ', id);
    end
    fprintf(fIds, '\n');
  end
  fclose(fIds);

  % % OBSOLETE FILES:
  % fNameMetacentrumTaskTemplate = [expDir filesep 'metacentrum_task_template.sh'];
  % fNameTask       = [experimentPath filesep exp_id '_task.sh'];
  % fNameTaskShell  = ['$EXPPATH_SHORT/$EXPID/${EXPID}_task.sh'];
  % fNameMng        = [experimentPath filesep exp_id '_manager.sh'];
  % fNameAllIds     = [experimentPath filesep 'allids.txt'];
  % copyfile(fNameMetacentrumTaskTemplate,        fNameTask);
  % if (isunix) fileattrib(fNameTask, '+x'); end
end
