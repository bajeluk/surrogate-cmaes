function expInitMetaLearn(exp_id)
% EXPINITMETALEARN -- Meta-learning experiment initialization
%
% EXP_ID -- string with an unique name of the experiment
%
% expInit('exp_name') has to be run before running the experiment.
% The experiment definitions are located in m-files 'exp/experiments/[exp_id].m'.
%
% SEE ALSO: generateStructOpts.m, getParamIndexVector.m, expUtilTest.m

  exppath_short = fileparts(mfilename('fullpath'));

  expFile = fullfile(exppath_short, 'experiments', [exp_id '.m']);
  if (~exist(expFile, 'file'))
    error('Error: experiment file "%s" does not exist.', expFile);
  end

  % load the core experiment configurations
  run(expFile);

%   % meta learning settings
%   metaParamDef = generateParamDefStruct(metaParams);

  % prepare a directory for experiments settings and results: exp/experiments/[exp_id]
  experimentPath = fullfile(exppath_short, 'experiments', exp_id);
  mkdir(experimentPath);
  
  % Model settings
  % -- enumerate all the possible combinations for each
  %    model type separately
  modelParamDef = generateParamDefStruct(modelOptions);

  for i = 1:length(modelParamDef)
    multiStruct = generateStructOpts(modelParamDef(i).values);
    modelParamDef(i).values = convertMultiStructToCell(multiStruct);
  end

  % save the definition of the experiment into 'scmaes_params.mat'
  save([experimentPath filesep 'metalearn_params.mat'], ...
    'modelParamDef', 'exp_id', 'exppath_short', 'exp_description', 'opts');

  % prepare other directories and files in the 'experimentPath'
%   params = [metaParamDef];
%   nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);

  % each model type has independent full factorial design
  get_name = @(x) x.name;
  get_val_len = @(x) length(x.values);
  modelNames = arrayfun(get_name, modelParamDef, 'UniformOutput', false);
  modelCombinations = arrayfun(get_val_len, modelParamDef);

  prepareMetacentrumFiles(exppath_short, experimentPath, exp_id, modelNames, ...
    modelCombinations);
end


function c = convertMultiStructToCell(st)
  c = cell(1,length(st));
  for i = 1:length(st)
    c{i} = st(i);
  end
end


function prepareMetacentrumFiles(expDir, experimentPath, exp_id, ...
  modelNames, modelCombinations)
  % expDir         = path to   'exp/'
  % experimentPath = path to   'exp/experiment/exp_id'
  % exp_id         = string with the unique experiment identifier

  % Export a textfile '<model>_ids.txt' with all the combinations' IDs for
  % each model
  for i = 1:length(modelNames)
    modelName = modelNames{i};
    n = modelCombinations(i);
    fNameModelIds = [experimentPath filesep modelName '_ids.txt'];
    fModelIds = fopen(fNameModelIds, 'w');
    if (fModelIds > 0)
      for id = 1:n, fprintf(fModelIds, '%d ', id); end
      fprintf(fModelIds, '\n');
    end
    fclose(fModelIds);
  end

  % Export a textfile 'allids.txt' with all the combinations' IDs
  fNameAllIds       = [experimentPath filesep 'allids.txt'];
  fIds = fopen(fNameAllIds, 'w');
  if (fIds > 0)
    for id = 1:sum(modelCombinations)
      fprintf(fIds, '%d ', id);
    end
    fprintf(fIds, '\n');
  end
  fclose(fIds);
end
