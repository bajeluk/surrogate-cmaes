function ds = modelTestSets(exp_id, fun, dim, inst, opts)
% ds = modelTestSets(exp_id, fun, dim) creates/loads dataset from 
% experiment.
%
% Input:
%   exp_id  - experiment id where the CMA-ES logs are saved
%             (there has to be only one variant!)
%   fun     - BBOB function numbers to load
%   dim     - dimensions to load (2, 3, 5, 10, ...)
%   maxEval - maximal number of evaluations times dimension to load
%
% Output:
%   ds - loaded data | #fun x #dim cell-array
%
% See Also:
%   datasetFromInstances, testModels, testOneModel

  if (~exist('opts', 'var'))
    opts = struct(); end
  if (~exist('exp_id', 'var'))
    exp_id = 'exp_doubleEC_21_log'; end
  if (~exist('fun', 'var'))
    fun = 1; end
  if (~exist('dim', 'var'))
    dim = 2; end
  if (~exist('inst', 'var'))
    inst = 1; end

  opts.outputDirname = defopts(opts, 'outputDirname', 'exp_GPtest_01');
  opts.datasetName   = defopts(opts, 'datasetName', 'DTS_005');
  opts.maxEval       = defopts(opts, 'maxEval', 250);
  % exp_id of the logging (DTS-)CMA-ES experiment
  opts.inputExp_id   = exp_id;
  % The number of generated datasets per one CMA-ES run
  opts.nSnapshotsPerRun = defopts(opts, 'nSnapshotsPerRun', 10);

  % set random seed due to reproducibility of default dataset
  rng(opts.maxEval)

  % path settings
  opts.exppath_short = fullfile('exp', 'experiments');
  opts.exppath = fullfile(opts.exppath_short, opts.inputExp_id);
  outputDir = fullfile('exp', 'experiments', opts.outputDirname);
  [~, ~] = mkdir(outputDir);
  outputDir = fullfile(outputDir, 'dataset');
  opts.datasetFile = fullfile(outputDir, [opts.datasetName '.mat']);

  % check experiment parameters
  opts.paramsMatFile = fullfile(opts.exppath, 'scmaes_params.mat');
  if exist(opts.paramsMatFile, 'file')
    exp_par = load(opts.paramsMatFile);
    exp_dim = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'dimensions')).values);
    exp_fun = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'functions')).values);
    exp_inst = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'instances')).values); 
  else
    warning('No scmaes_params.mat found in %s directory. Using default fun and dim settings.', opts.exppath)
    exp_dim = 2;
    exp_fun = 1;
    exp_inst = 1;
  end

  % take only dimensions and function which exist in the logging experiment
  dim = restricToDataset(dim, exp_dim, 'Dimensions');
  fun = restricToDataset(fun, exp_fun, 'Functions');
  inst = restricToDataset(inst, exp_inst, 'Instances');

  assert(~isempty(inst), 'There are no instances to process. Exitting.');
  
  if (exist(opts.datasetFile, 'file'))
    fprintf('The dataset file "%s" already existed.\n   Copying to "%s.bak".\n', opts.datasetFile, opts.datasetFile);
    copyfile(opts.datasetFile, [opts.datasetFile '.bak']);
    f_ds = load(opts.datasetFile);
    if (isfield(f_ds, 'ds') && ~isempty(f_ds.ds))
      ds = f_ds.ds;
      dim = f_ds.dim;
      fun = f_ds.fun;
      inst = f_ds.inst;
    end
  else
    f_ds = struct('ds', {{}}, 'fun', {[]}, 'dim', {[]}, 'inst', {[]}, 'maxEval', {opts.maxEval});
    ds = cell(length(fun), length(dim), length(inst));
  end

  % dimension loop
  for di = 1:length(dim)
    d_exp = find(dim(di) == exp_dim);

    % function loop
    for fi = 1:length(fun)
      f_exp = find(fun(fi) == exp_fun);
      id = (d_exp-1)*length(exp_fun) + f_exp;
      fprintf('#### f%d in %dD ####\n', fun(fi), dim(di));

      % instances loop
      instancesDone = false(1, length(inst));
      for ii = 1:length(inst)
        i_exp = find(inst(ii) == exp_inst);
        if (~isempty(ds{fi, di, ii}) && isstruct(ds{fi, di, ii}))
          instancesDone(ii) = true;
        end
      end
      if (all(instancesDone))
        fprintf('All instances are done. Advancing to the next function...\n');
        continue
      end

      % load dataset from saved modellog/cmaes_out of the corresponding instance
      ds_actual = datasetFromInstances(opts, opts.nSnapshotsPerRun, ...
          fun(fi), dim(di), inst, id);
        
      ds(f_exp, d_exp, :) = ds_actual;

    % function loop end  
    end
  % dimension loop end
  end

  % save default dataset
  fun = union(fun, f_ds.fun);
  dim = union(dim, f_ds.dim);
  inst = union(inst, f_ds.inst);
  maxEval = opts.maxEval;
  save(opts.datasetFile, 'ds', 'fun', 'dim', 'inst', 'maxEval')
end

function [intToTest] = restricToDataset(intToTest, dataInt, identifier)
  id_ok = ismember(intToTest, dataInt);
  if ~all(id_ok)
    warning('%s %s are not in the experiment', identifier, ...
      strjoin(arrayfun(@num2str, intToTest(~id_ok), 'UniformOutput', false), ','))
    intToTest = intToTest(id_ok);
  end
end
