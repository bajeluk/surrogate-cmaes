function status = metacentrum_testmodels(exp_id, exppath_short, func_str, dim_str, inst_str, opts_str, dataset)
% exp_GPtest_01 model testing experiment -- Matlab part to be MCR-compiled
%
% Usage:
%   metacentrum_testmodels(exp_id, exppath_short, func_str, dim_str, inst_str, opts_str, dataset)

  % Input parameter settings
  %
  % FUN/DIM/INST input parse settings
  if (~exist('func_str', 'var'))
    func_str = []; end
  if (~exist('dim_str', 'var'))
    dim_str = []; end
  if (~exist('inst_str', 'var'))
    inst_str = []; end
  if (~exist('opts_str', 'var'))
    opts_str = []; end

  func          = parseCmdParam('func_str', func_str, 1:24);
  dims          = parseCmdParam('dim_str',  dim_str,  [2, 5, 10]);
  instances     = parseCmdParam('inst_str', inst_str, [1:5, 41:50]);
  cmd_opts      = parseCmdParam('opts_str', opts_str, struct());

  % run the script EXP_ID.m if it exists
  opts = struct();
  expScript = fullfile(exppath_short, [exp_id '.m']);
  if (exist(expScript, 'file'))
    run(expScript);
  end
  % and re-write the options specified on command-line
  if (isstruct(cmd_opts) && ~isempty(cmd_opts))
    cmd_fnames = fieldnames(cmd_opts);
    for i = 1:length(cmd_fnames)
      opts.(cmd_fnames{i}) = cmd_opts.(cmd_fnames{i});
    end
  end

  % EXPID -- unique experiment identifier (directory with this name
  %          is expected to hold the dataset and the results will be
  %          placed there, too)
  opts.exp_id     = exp_id;
  % EXPPATH_SHORT
  opts.exppath_short = exppath_short;

  % dataset name (filename w/o extension of the dataset) | string
  % or struct with the field '.ds' with 3D cell array with the data for {d, f, i}'s
  if (~exist('dataset', 'var') || isempty(dataset))
    opts.dataset  = defopts(opts, 'dataset',  'DTS_005');
  else
    opts.dataset = dataset;
  end

  % type of model to test according to the ModelFactory | string
  opts.modelType = defopts(opts, 'modelType', 'gp');

  % statistics to compute
  opts.statistics = defopts(opts, 'statistics', ...
      { 'mse', 'mzoe', 'kendall', 'rankmse', 'rankmzoe', 'rde' });

  % Maximal number of function evaluation per dimension to consider
  opts.maxEvals = defopts(opts, 'maxEvals', 250);

  % path settings
  opts.scratch = defopts(opts, 'scratch', getenv('SCRATCHDIR'));

  % other settings
  opts.rewrite_results = defopts(opts, 'rewrite_results', false);

  % directory with the dataset and results
  opts.exppath = fullfile(opts.exppath_short, opts.exp_id);
  % specifying the dataset -- expand the filename if dataset is string and not file
  if (ischar(opts.dataset) && ~exist(opts.dataset, 'file'))
    opts.dataset = fullfile(opts.exppath, 'dataset', [opts.dataset, '.mat']);
  end
  % defSetFile = fullfile(scratch, 'model', 'defData', ['defSet_', num2str(maxEvals), 'FE.mat']);

  % use this default model options if not specified in the experiment
  % script exp/experiments/EXP_ID.m
  if (~exist('modelOptions', 'var'))
    modelOptions.useShift        = false;
    modelOptions.predictionType  = 'sd2';
    modelOptions.trainAlgorithm  = 'fmincon';
    modelOptions.covFcn          = '{@covMaterniso, 5}';
    modelOptions.normalizeY      = true;
    modelOptions.hyp.lik         = log(0.01);
    modelOptions.hyp.cov         = log([0.5; 2]);
    modelOptions.covBounds       = [ [-2;-2], [25;25] ];
    modelOptions.likBounds       = log([1e-6, 10]);

    % Full factorial design of the following parameters
    modelOptions.trainsetType    = { 'nearestToPopulation', 'nearest', 'clustering', 'allPoints' };
    modelOptions.trainRange      = { 1.0, 0.999 };
    modelOptions.trainsetSizeMax = { '5*dim', '10*dim', '15*dim', '20*dim' };
    modelOptions.meanFcn         = { 'meanConst', 'meanLinear' };
    modelOptions.covFcn          = { '{@covSEiso}', '{@covSEard}', ...
                                '{@covMaterniso, 5}', '{@covMaterniso, 3}' };
  end

  % modelOptions can be specified in opts_str, EXP_ID.m, or the above
  % defined structure will be used
  modelOptions = defopts(opts, 'modelOptions', modelOptions);

  % Combine options into full factorial design
  modelOptions_fullfact  = combineFieldValues(modelOptions);

  % restrict the fullfactorial design only to some indices
  % if specified in opts
  if (isfield(opts, 'modelOptionsIndices') && ~isempty(opts.modelOptionsIndices))
    opts.modelOptionsIndices = myeval(opts.modelOptionsIndices);
    modelOptions_fullfact = modelOptions_fullfact(opts.modelOptionsIndices);
  end

  %% create testing dataset
  %ds = modelTestSets('exp_doubleEC_21_log15', func, dims, instances, opts);

  fprintf('== Summary of the testing assignment ==\n');
  fprintf('   # of models:  %d\n', length(modelOptions_fullfact));
  fprintf('   functions:    %s\n', num2str(func));
  fprintf('   dimensions:   %s\n', num2str(dims));
  fprintf('   instances:    %s\n', num2str(instances));
  fprintf('=======================================\n');

  %% test chosen models
  modelFolders = testModels(modelOptions_fullfact, opts, func, dims, instances);

  %% load and calculate results
  [rdeTable, mseTable, RDEs, MSEs] = modelStatistics(modelFolders, func, dims, instances);

  save(fullfile(opts.exppath, 'stats.mat'), 'modelFolders', 'rdeTable', 'mseTable', 'RDEs', 'MSEs');

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
