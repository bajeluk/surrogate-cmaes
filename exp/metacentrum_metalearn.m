function status = metacentrum_metalearn(exp_id, exppath_short, dim_str, func_str, ...
  inst_str, size_str, design_str, model_str, opts_str, dataset_path)
%METACENTRUM_METALEARN Top-level function for running meta learning experiments. Supposed to be
%   run as a MCR-compiled binary.
%   METACENTRUM_METALEARN(exp_id, exppath_short, func_str, dim_str, inst_str, opts_str, dataset_path)
%      exp_id        -- experiment id
%      exppath_short -- relative path to experiment folder
%      dim_str       -- dimensinalities to test (string)
%      func_str      -- function ids list (string)
%      inst_str      -- instances to test (string)
%      size_str      -- dataset sizes to test (string)
%      design_str    -- design types to test (string)
%      model_str     -- model types to test (string)
%      opts_str      -- options with model parameters grid
%      dataset_path  -- path to metalearning data sets

  % Input parameter settings
  %
  % FUN/DIM/INST/... input parse settings
  if (~exist('func_str', 'var')), func_str = []; end
  if (~exist('dim_str', 'var')), dim_str = []; end
  if (~exist('inst_str', 'var')), inst_str = []; end
  if (~exist('size_str', 'var')), size_str = []; end
  if (~exist('design_str', 'var')), design_str = []; end
  if (~exist('model_str', 'var')), model_str = []; end
  if (~exist('opts_str', 'var')), opts_str = []; end

  dims          = parseCmdParam('dim_str',  dim_str,  [2, 5, 10, 20]);
  func          = parseCmdParam('func_str', func_str, 1:24);
  instances     = parseCmdParam('inst_str', inst_str, [1:5, 41:50]);
  Ns            = parseCmdParam('size_str', size_str, {'50 * dim'});
  designs       = parseCmdParam('design_str', design_str, {'lhs'});
  models        = parseCmdParam('model_str', model_str, {'gp', 'rf'});
  cmd_opts      = parseCmdParam('opts_str', opts_str, struct());

  % load options and full factorial designs for model settings
  params_file = fullfile(exppath_short, exp_id, 'metalearn_params.mat');
  load(params_file, 'opts', 'modelParamDef');

  % and re-write the options specified on command-line
  if (isstruct(cmd_opts) && ~isempty(cmd_opts))
    cmd_fnames = fieldnames(cmd_opts);
    for i = 1:length(cmd_fnames)
      opts.(cmd_fnames{i}) = cmd_opts.(cmd_fnames{i});
    end
  end

  % EXPID -- unique experiment identifier
  % results will be placed into this dir
  opts.exp_id     = exp_id;
  % EXPPATH_SHORT
  opts.exppath_short = exppath_short;

  % dataset path name (filename w/o extension of the dataset) | string
  if (~exist('dataset_path', 'var') || isempty(dataset_path))
    opts.dataset_path = defopts(opts, 'dataset_path',  'data_metalearning');
  else
    opts.dataset_path = dataset_path;
  end

  % path settings
  opts.scratch = defopts(opts, 'scratch', getenv('SCRATCHDIR'));

  % other settings
  opts.rewrite_results = defopts(opts, 'rewrite_results', false);

  % directory with the dataset and results
  opts.exppath = fullfile(opts.exppath_short, opts.exp_id);

  % specifying the dataset -- expand the filename if dataset is string and not file
  if (ischar(opts.dataset_path) && ~exist(opts.dataset_path, 'dir'))
    opts.dataset_path = fullfile(opts.scratch, opts.dataset_path);
  end

  % Load a full factorial design for each model
  modelOptions_fullfact = struct();
  
  for i = 1:length(modelParamDef)
    name = modelParamDef(i).name;
    values = modelParamDef(i).values;
    modelOptions_fullfact.(name) = values;
  end

  nOptions = 0;
  for model_cell = models
    model = model_cell{:};
    if ~isfield(modelOptions_fullfact, model)
      error('Options for model ''%s'' not defined', model);
    end
    nOptions = nOptions + length(modelOptions_fullfact.(model));
  end
  
%   % restrict the full factorial design only to some indices
%   % if specified in opts
%   if (isfield(opts, 'modelOptionsIndices') && ~isempty(opts.modelOptionsIndices))
%     opts.modelOptionsIndices = myeval(opts.modelOptionsIndices);
%     modelOptions_fullfact = modelOptions_fullfact(opts.modelOptionsIndices);
%   end

  %% create testing dataset
  fprintf('== Summary of the testing assignment ==\n');
  fprintf('   # of models:  %d\n', nOptions);
  fprintf('   functions:    %s\n', num2str(func));
  fprintf('   dimensions:   %s\n', num2str(dims));
  fprintf('   instances:    %s\n', num2str(instances));
  fprintf('   sizes:        %s\n', strjoin(Ns));
  fprintf('   designs:      %s\n', strjoin(designs));
  fprintf('   models:       %s\n', strjoin(models));
  fprintf('=======================================\n');

  %% test chosen models
  %modelFolders = testModels(modelOptions_fullfact, opts, func, dims, instances);

  status = 0;
  return;
end

function out = parseCmdParam(name, value, defaultValue)
  if (isempty(value))
    out = defaultValue;
  elseif (ischar(value))
    out = myeval(value);
  else
    error('%s has to be string for eval()-uation', name);
  end
end
