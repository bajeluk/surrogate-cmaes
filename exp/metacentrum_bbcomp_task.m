function status = metacentrum_bbcomp_task(exp_id, exppath_short, problemID_str, opts_str)
% metacentrum_bbcomp_task -- Matlab part of BBComp 2017 competition scripts to be MCR-compiled
%
% Usage:
%   metacentrum_bbcomp_task(exp_id, exppath_short, problemID_str, opts_str)

  % Input parameter settings
  %
  % ID/opts input parse settings
  if (~exist('problemID_str', 'var'))
    problemID_str = []; end
  if (~exist('opts_str', 'var'))
    opts_str = []; end

  id            = parseCmdParam('problemID_str', problemID_str, 1:5);
  cmd_opts      = parseCmdParam('opts_str', opts_str, struct());

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
  if (isstruct(cmd_opts) && ~isempty(cmd_opts))
    cmd_fnames = fieldnames(cmd_opts);
    for i = 1:length(cmd_fnames)
      opts.(cmd_fnames{i}) = cmd_opts.(cmd_fnames{i});
    end
  end

  % EXPID -- unique experiment identifier (directory in where 
  %          the results are expected to be placed)
  opts.exp_id     = exp_id;
  % EXPPATH_SHORT
  opts.exppath_short = exppath_short;
  opts.exppath = [exppath_short filesep exp_id];
  RESULTSFILE = [opts.exppath '/' exp_id '_results_' num2str(id) '.mat'];
  if (isempty(getenv('SCRATCHDIR')))
    OUTPUTDIR = '/tmp/job_output';     % set OUTPUTDIR empty if $SCRATCHDIR var does not exist
  else
    OUTPUTDIR = [getenv('SCRATCHDIR') filesep 'job_output'];
  end

  % datapath is typically on the computing node in   $SCRATCHDIR/job_output/bbcomp_output:
  datapath = [OUTPUTDIR filesep 'bbcomp_output'];
  if (isempty(strfind(datapath, exppath_short)))
    % localDatapath is typically on the NFS server in
    %   $HOME/prg/surrogate-cmaes/exp/experiments/$EXPID/bbcomp_output_tmp
    % and the bbcomp results are stored here after each completed instance
    localDatapath = [opts.exppath filesep 'bbcomp_output_tmp'];
    [~, ~] = mkdir(localDatapath);
  end
  [~, ~] = mkdir(datapath);
  surrogateParams.datapath = datapath;

  % Metacentrum task and node properties
  cmd_opts.machine = '';
  nodeFile = fopen(getenv('PBS_NODEFILE'), 'r');
  if (nodeFile > 0)
    cmd_opts.machine = fgetl(nodeFile);
    fclose(nodeFile);
  end

  fprintf('== Summary of the testing assignment ==\n');
  fprintf('   problemID:    %s\n', num2str(id));
  fprintf('=======================================\n');
  % Matlab should have been called from a SCRACHDIR
  startup;

  t0 = clock;
  % Initialize random number generator
  exp_settings.seed = myeval(defopts(cmd_opts, 'seed', 'floor(sum(100 * clock()))'));
  rng(exp_settings.seed);

  %
  %
  % the computation itself
  %

  % MOCK IMPLEMENTATION
  FUN = @frosen;
  dim = 2;
  maxfunevals = 5000;

  surrogateParams.archive = Archive(dim);
  surrogateParams.startTime = tic;
  nEvals = 0;
  xstart = rand(dim, 1);

  while nEvals < maxfunevals

    [x, y_evals, stopflag, archive, varargout] = opt_s_cmaes_bbcomp(FUN, dim, ...
        maxfunevals, cmaesParams, surrogateParams, xstart);

    nEvals = size(archive.X, 1);

    xstart = cmaesRestartPoint(archive);
  end


  % copy the bbcomp results onto persistant storage if outside EXPPATH
  if (~isempty(OUTPUTDIR) && ~strcmpi(OUTPUTDIR, opts.exppath) && isunix)
    % copy the output to the final storage (if OUTPUTDIR and EXPPATH differs)
    system(['cp -pR ' OUTPUTDIR '/* ' opts.exppath '/']);
  end

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
