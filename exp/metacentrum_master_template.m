function job = metacentrum_master_template(exp_id, varargin)
% job = metacentrum_master_template(exp_id, task_id, walltime) runs chosen
% tasks of experiment exp_id on metacentrum.
%
% Input:
%   exp_id   - experiment ID string
%   task_id  - vector of integer ID's of the parameter-combinations to try 
%            - default: ALL of the combinations
%   walltime - string defining maximum (wall)time for Metacentrum machines
%            - default: '4h'
%
% See Also:
%   metacentrum_task_matlab

  if nargin == 0
    help metacentrum_master_template
    if nargout > 0
      job = [];
    end
    return
  end
  
  pathstr = [fileparts(mfilename('fullpath')) filesep 'experiments'];
  exppath = [pathstr filesep exp_id];
  load([exppath filesep 'scmaes_params.mat']);
  % this is old version with a fixed path :(
  % load(['/storage/plzen1/home/' getenv('LOGNAME') '/prg/surrogate-cmaes/exp/experiments/' exp_id '/scmaes_params.mat']);

  fNameJobDef = [exppath_short filesep exp_id filesep exp_id '_metajob.mat'];
  cd([exppath_short filesep '..' filesep '..']);
  % BE CAREFUL!! this is manual PATH addition
  addpath('exp/util');  % for getParamsFromIndex()
  addpath('src/util');  % for structReduce()
  addpath([exppath_short filesep '..']);        % for metacentrum_task_matlab()
  if (nargin >= 2 && ~isempty(varargin{1}))
    combsToRun = varargin{1};
  else
    params = [bbParamDef, sgParamDef, cmParamDef];
    nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);
    combsToRun = 1:nCombinations;
  end
  if (nargin >= 3 && ~isempty(varargin{2}))
    walltime = varargin{2};
  else
    walltime = '4h';
  end

  pbs_max_workers = 50;
  pbs_params = ['-l walltime=' walltime ',nodes=^N^:ppn=1,mem=1gb,scratch=1gb,matlab_MATLAB_Distrib_Comp_Engine=^N^'];

  % sched = findResource('scheduler','type','torque');
  % set(sched,'ClusterOsType', 'unix');
  % set(sched,'DataLocation',['/storage/plzen1/home/' getenv('LOGNAME') '/matlab']);
  % set(sched,'HasSharedFilesystem', true);
  % set(sched,'ClusterMatlabRoot',matlabroot); % v celem MetaCentru je MATLAB dostupny ze stejneho adresare
  % % set(sched,'SubmitArguments','-q short');
  % set(sched,'ResourceTemplate', pbs_params);
  % get(sched)

  while 1
    [tf msg] = license('checkout','Distrib_Computing_Toolbox');
    if tf==1, break, end
    display(strcat(datestr(now),' waiting for licence '));
    pause(4);
  end

  cl = parallel.cluster.Torque;
  pause(2);
  [~, ~] = mkdir(exppath_short, '../matlab_jobs')
  cl.JobStorageLocation = [exppath_short filesep '../matlab_jobs'];
  cl.ClusterMatlabRoot = matlabroot;
  cl.OperatingSystem = 'unix';
  cl.ResourceTemplate = pbs_params;
  cl.HasSharedFilesystem = true;
  cl.NumWorkers = pbs_max_workers;

  job = createJob(cl);

  i = 1;
  for id = combsToRun
    [bbParams, sgParams] = getParamsFromIndex(id, bbParamDef, sgParamDef, cmParamDef);
    fileID = [num2str(bbParams.functions(end)) '_' num2str(bbParams.dimensions(end)) 'D_' num2str(id)];
    if (isfield(sgParams, 'modelType'))
      model = sgParams.modelType;
    else
      model = 'NONE';
    end

    metaOpts.logdir = logDir;
    metaOpts.model = model;
    metaOpts.nInstances = length(bbParams.instances);
    fprintf('Setting up job ID %d / %d (f%d/%dD)...\n', id, length(combsToRun), bbParams.functions(end), bbParams.dimensions(end));
    tasks(i) = createTask(job, @metacentrum_task_matlab, 0, {exp_id, exppath_short, id, metaOpts});
    i = i + 1;
  end

  tasks

  submit(job)

  save(fNameJobDef, 'exp_id', 'metaOpts', 'combsToRun', 'pbs_params');
end
