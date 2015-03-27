function metacentrum_master_template(exp_id)
  load(['/storage/plzen1/home/' getenv('LOGNAME') '/prg/surrogate-cmaes/exp/experiments/' exp_id '/scmaes_params.mat']);

  fNameJobDef = [exppath_short filesep exp_id filesep exp_id '_metajob.mat'];
  cd([exppath_short filesep '..' filesep '..']);
  % BE CAREFUL!! this is manual PATH addition
  addpath('exp/util');  % for getParamsFromIndex()
  addpath('src/util');  % for structReduce()
  addpath([exppath_short filesep '..']);        % for metacentrum_task_matlab()
  params = [bbParamDef, sgParamDef, cmParamDef];
  nCombinations = structReduce(params, @(s,x) s*length(x.values), 1);

  pbs_max_workers = 30;
  pbs_params = '-l walltime=4h,nodes=^N^:ppn=1,mem=2gb,scratch=1gb,matlab=1,matlab_Statistics_Toolbox=1,matlab_Optimization_Toolbox=1,matlab_MATLAB_Distrib_Comp_Engine=^N^';

  % sched = findResource('scheduler','type','torque');
  % set(sched,'ClusterOsType', 'unix');
  % set(sched,'DataLocation',['/storage/plzen1/home/' getenv('LOGNAME') '/matlab']);
  % set(sched,'HasSharedFilesystem', true);
  % set(sched,'ClusterMatlabRoot',matlabroot); % v celem MetaCentru je MATLAB dostupny ze stejneho adresare
  % % set(sched,'SubmitArguments','-q short');
  % set(sched,'ResourceTemplate', pbs_params);
  % get(sched)

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

  for id = 1:nCombinations
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
    fprintf('Setting up job ID %d (f%d/%dD)...\n', id, bbParams.functions(end), bbParams.dimensions(end));
    tasks(id) = createTask(job, @metacentrum_task_matlab, 0, {exp_id, exppath_short, bbParams.functions(end), bbParams.dimensions(end), id, metaOpts});
  end

  tasks

  submit(job)

  save(fNameJobDef);
end
