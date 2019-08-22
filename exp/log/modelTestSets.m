function ds = modelTestSets(exp_id, fun, dim, inst, opts)
% ds = modelTestSets(exp_id, fun, dim, inst, opts) creates/loads dataset
% from experiment 'exp_id'.
%
% Input:
%   exp_id  - experiment id where the CMA-ES logs are saved
%   fun     - BBOB function numbers to load
%   dim     - dimensions to load (2, 3, 5, 10, ...)
%   inst    - instances to load (1:5, ...)
%   opts    - struct with additional settings as fields:
%     'datasetName'          - name of resulting dataset | 'DTS_005'
%     'isForModelPool'       - dataset for ModelPool testing | false
%     'loadModels'           - load saved models for ModelPool | false
%     'maxEval'              - maximal number of evaluations times
%                              dimension to load | 250
%     'nPreviousGenerations' - number of previous generations to load for
%                              ModelPool testing | 0
%     'nSnapshotsPerRun'     - number of generated datasets per one CMA-ES
%                              run | 10
%     'outputDirname'        - name of dataset output directory |
%                              'exp_modeltest_01'
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

  opts.outputDirname = defopts(opts, 'outputDirname', 'exp_modeltest_01');
  opts.datasetName   = defopts(opts, 'datasetName', 'DTS_005');
  opts.maxEval       = defopts(opts, 'maxEval', 250);
  % exp_id of the logging (DTS-)CMA-ES experiment
  opts.inputExp_id   = exp_id;
  % The number of generated datasets per one CMA-ES run
  opts.nSnapshotsPerRun = defopts(opts, 'nSnapshotsPerRun', 10);
  % dataset for ModelPool testing needs additional data
  opts.isForModelPool = defopts(opts, 'isForModelPool', false);
  opts.nPreviousGenerations = defopts(opts, 'nPreviousGenerations', 0);
  opts.loadModels     = defopts(opts, 'loadModels', false);

  % set random seed due to reproducibility of default dataset
  rng(opts.maxEval)

  % path settings
  opts.exppath_short = fullfile('exp', 'experiments');
  opts.exppath = fullfile(opts.exppath_short, opts.inputExp_id);
  outputDir = fullfile('exp', 'experiments', opts.outputDirname);
  [~, ~] = mkdir(outputDir);
  outputDir = fullfile(outputDir, 'dataset');
  [~, ~] = mkdir(outputDir);
  opts.datasetFile = fullfile(outputDir, [opts.datasetName '.mat']);

  % check experiment parameters
  opts.paramsMatFile = fullfile(opts.exppath, 'scmaes_params.mat');
  if exist(opts.paramsMatFile, 'file')
    exp_par = load(opts.paramsMatFile);
    exp_dim = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'dimensions')).values);
    exp_fun = cell2mat(exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'functions')).values);
    exp_inst_cell = exp_par.bbParamDef(strcmp({exp_par.bbParamDef.name}, 'instances')).values;
    exp_inst = cell2mat(exp_inst_cell);
    exp_model_settings = 1:prod(arrayfun(@(s) length(s.values), exp_par.sgParamDef));
  else
    warning('No scmaes_params.mat found in %s directory. Using default fun and dim settings.', opts.exppath)
    exp_dim = 2;
    exp_fun = 1;
    exp_inst_cell = { 1 };
    exp_inst = 1;
    exp_model_settings = 1;
  end

  % TODO: make this configurable
  modelSettings = exp_model_settings;

  % take only dimensions and function which exist in the logging experiment
  dim = restricToDataset(dim, exp_dim, 'Dimensions');
  fun = restricToDataset(fun, exp_fun, 'Functions');
  inst = restricToDataset(inst, exp_inst, 'Instances');
  % TODO: should also be user-specified and restricted to dataset
  modelSettings = restricToDataset(modelSettings, exp_model_settings, 'Model settings');

  IDsPerInstance = ...
    (prod(arrayfun(@(s) length(s.values), exp_par.bbParamDef)) / (length(exp_fun)*length(exp_dim))) * ...
    prod(arrayfun(@(s) length(s.values), exp_par.cmParamDef)) * ...
    prod(arrayfun(@(s) length(s.values), exp_par.sgParamDef));

  assert(~isempty(inst), 'There are no instances to process. Exitting.');

  fprintf('== Summary of the dataset being created ==\n');
  fprintf('   functions:    %s\n', num2str(fun));
  fprintf('   dimensions:   %s\n', num2str(dim));
  fprintf('   instances:    %s\n', num2str(inst));
  fprintf('==========================================\n');

  is_ds_loaded = false;
  if (exist(opts.datasetFile, 'file'))
    fprintf('The dataset file "%s" already existed.\n   Copying to "%s.bak".\n', opts.datasetFile, opts.datasetFile);
    copyfile(opts.datasetFile, [opts.datasetFile '.bak']);
    f_ds = load(opts.datasetFile);
    if (isfield(f_ds, 'ds') && ~isempty(f_ds.ds))
      is_ds_loaded = true;
    end
  else
    f_ds = struct('ds', {{}}, 'fun', {[]}, 'dim', {[]}, 'inst', {[]}, 'modelSettings', {[]}, 'maxEval', {opts.maxEval});
  end

  % prepare output dataset
  ds = cell(length(fun), length(dim), length(inst), length(exp_model_settings));

  % dimension loop
  for di = 1:length(dim)
    d_exp = find(dim(di) == exp_dim, 1);

    % function loop
    for fi = 1:length(fun)
      f_exp = find(fun(fi) == exp_fun, 1);

      % check whether if we have some instances loaded
      instancesDone = false(1, length(inst));
      if (is_ds_loaded)
        for ii = 1:length(inst)
          % i_exp = find(inst(ii) == exp_inst, 1);
          d_loaded = find(dim(di) == f_ds.dim, 1);
          f_loaded = find(fun(fi) == f_ds.fun, 1);
          i_loaded = find(inst(ii) == f_ds.inst, 1);

          if ((~isempty(d_loaded) && ~isempty(f_loaded) && ~isempty(i_loaded)) ...
              && ~any(cellfun(@isempty, f_ds.ds(f_loaded, d_loaded, i_loaded, exp_model_settings))) ...
              && all(cellfun(@isstruct, f_ds.ds(f_loaded, d_loaded, i_loaded, exp_model_settings))))
            instancesDone(ii) = true;

            ds(fi, di, ii, exp_model_settings) = f_ds.ds(f_loaded, d_loaded, i_loaded, exp_model_settings);
          end
        end
        if (all(all(instancesDone)))
          fprintf('All instances are done. Advancing to the next function...\n');
          continue
        end
      end

      % instances loop
      for ii_cell = 1:length(exp_inst_cell)
        id = (d_exp-1)*length(exp_fun)*length(exp_inst_cell)*IDsPerInstance + ...
          (f_exp-1)*length(exp_inst_cell)*IDsPerInstance + ...
          (ii_cell-1)*length(exp_inst_cell) + ...
          exp_model_settings;
        fprintf('#### f%d in %dD (id=%s) ####\n', fun(fi), dim(di), num2str(id));

        instancesInThisID = exp_inst_cell{ii_cell};

        % load dataset from saved modellog/cmaes_out of the corresponding instance
        [~, instanceIndicesToProcess] = ismember(instancesInThisID, inst);
        instanceIndicesToProcess(instanceIndicesToProcess == 0) = [];
        instancesToProcess = intersect(inst(instanceIndicesToProcess), inst(~instancesDone));
        [~, instanceIndicesToProcess] = ismember(instancesToProcess, inst);

        if (~isempty(instancesToProcess))
          % create #inst x #id dataset (cell-array of structures with
          % DTS-CMA-ES state variables as fields)
          ds_actual = datasetFromInstances(exp_id, fun(fi), dim(di), inst(instanceIndicesToProcess), id, opts);
          ds(fi, di, instanceIndicesToProcess, exp_model_settings) = ds_actual;
        end
      end  % instances loop end
    end  % function loop end

    % save the dataset
    maxEval = opts.maxEval;
    nSnapshotsPerRun = opts.nSnapshotsPerRun;
    save(opts.datasetFile, 'ds', 'fun', 'dim', 'inst', 'modelSettings', 'maxEval', 'nSnapshotsPerRun', 'opts');

  end  % dimension loop end
end

function [intToTest] = restricToDataset(intToTest, dataInt, identifier)
  id_ok = ismember(intToTest, dataInt);
  if ~all(id_ok)
    warning('%s %s are not in the experiment', identifier, ...
      strjoin(arrayfun(@num2str, intToTest(~id_ok), 'UniformOutput', false), ','))
    intToTest = intToTest(id_ok);
  end
end
