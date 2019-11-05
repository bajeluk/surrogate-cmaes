function status = metacentrum_modeltestsets(exp_id, exppath_short, func_str, dim_str, inst_str, exp_log, dataset)
% metacentrum_modeltestsets - metacentrum Matlab part of modelTestSets 
% to be MCR-compiled.
%
% status = metacentrum_modeltestsets(exp_id, exppath_short, func_str, 
%   inst_str, ids_str, opts_str, exp_log, dataset)
%
% Input:
%   exp_id - unique experiment identifier | string
%   exppath_short - path to location of experiment folder | string
%   func_str      - string of functions to calculate | string
%   dim_str       - string of dimensions to calculate | string
%   inst_str      - string of instances to calculate | string
%   ids_str       - string of ids to calculate | string
%   opts_str      - string of options to calculate | string
%   exp_log       - logging experiment name (data source) | string
%   dataset       - dataset name | string
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
  % if (~exist('opts_str', 'var'))
  %   opts_str = []; end
  % if (~exist('ids_str', 'var'))
  %   ids_str = []; end

  func          = parseCmdParam('func_str', func_str, 1:24);
  dims          = parseCmdParam('dim_str',  dim_str,  [2, 5, 10]);
  instances     = parseCmdParam('inst_str', inst_str, [1:5, 41:50]);
  % cmd_opts      = parseCmdParam('opts_str', opts_str, struct());
  % ids           = parseCmdParam('ids_str', ids_str, 1);

  % run the script EXP_ID.m if it exists
  opts = struct();
  expScript = fullfile(exppath_short, [exp_id '.m']);
  if (exist(expScript, 'file'))
    % eval the script line by line (as it cannot be run()-ed when deployed)
    fid = fopen(expScript);
    tline = fgetl(fid);
    while (ischar(tline))
      eval(tline);
      tline = fgetl(fid);
    end
    fclose(fid);
  end
  % and re-write the options specified on command-line
  % if (isstruct(cmd_opts) && ~isempty(cmd_opts))
  %   cmd_fnames = fieldnames(cmd_opts);
  %   for i = 1:length(cmd_fnames)
  %     opts.(cmd_fnames{i}) = cmd_opts.(cmd_fnames{i});
  %   end
  % end

  % EXPID -- unique experiment identifier (directory with this name
  %          is expected to hold the results)
  opts.exp_id     = exp_id;
  % EXP_LOG -- logging experiment identifier (directory with this name
  %            is expected to hold the results of logging experiments)
  if (~exist('exp_log', 'var') || isempty(exp_log))
    exp_log = defopts(opts, 'exp_log', exp_id);
  end
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
  opts.mfts_only = defopts(opts, 'mfts_only', false);

  % directory with the dataset and results
  opts.exppath = fullfile(opts.exppath_short, opts.exp_id);
  % specifying the dataset -- expand the filename if dataset is string and not file
  if (ischar(opts.dataset) && ~exist(opts.dataset, 'file'))
    opts.dataset = fullfile(opts.exppath, 'dataset', [opts.dataset, '.mat']);
  end

  % export summary
  fprintf('== Summary of the testing assignment ==\n');
  fprintf('   functions:    %s\n', num2str(func));
  fprintf('   dimensions:   %s\n', num2str(dims));
  fprintf('   instances:    %s\n', num2str(instances));
  fprintf('=======================================\n');

  % create model test sets
  modelTestSets(exp_log, func, dims, instances, opts);

  % return status zero
  status = 0;
  return;
end

function out = parseCmdParam(name, value, defaultValue)
  if (isempty(value))
    out = defaultValue;
  elseif (ischar(value))
    fprintf('Run: parseCmdParam: value = %s\n', value)
    out = myeval(value);
  else
    error('%s has to be string for eval()-uation', name);
  end
end
