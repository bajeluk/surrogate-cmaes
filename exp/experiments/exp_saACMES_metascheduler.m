function exp_saACMES_metascheduler(varargin)
% single purpose function for setting up Ilya Loshchilov's saACM-ES test
% on Metacentrum
% splitting the tasks is done per function basis
% it is first run for dimension 20, and afterwords for dimensions [2,3,5,10]
% OPTIONAL PARAMETER: Metacentrum queue (2h, 4h, 1d, 2d, 4d...)

  functions = [2:24];
  dimensions = [20];

  disp('Current path of the script:');
  cwd = fileparts(mfilename('fullpath'))
  cd(cwd);
  % BE CAREFUL!! this is manual PATH addition
  addpath('saACMESlambdaRevMinIter3v2');  % for getParamsFromIndex()
  if (nargin >= 1 && ~isempty(varargin{1}))
    walltime = varargin{1};
  else
    walltime = '1d';
  end

  pbs_max_workers = 50;
  pbs_params = ['-l walltime=' walltime ',nodes=^N^:ppn=1,mem=1gb,scratch=1gb,matlab_MATLAB_Distrib_Comp_Engine=^N^'];

  while 1
    [tf msg] = license('checkout','Distrib_Computing_Toolbox');
    if tf==1, break, end
    display(strcat(datestr(now),' waiting for licence '));
    pause(4);
  end

  cl = parallel.cluster.Torque;
  pause(2);
  matlabjobsdir = [cwd '/../matlab_jobs'];
  [~, ~] = mkdir(matlabjobsdir)
  cl.JobStorageLocation = matlabjobsdir;
  cl.ClusterMatlabRoot = matlabroot;
  cl.OperatingSystem = 'unix';
  cl.ResourceTemplate = pbs_params;
  cl.HasSharedFilesystem = true;
  cl.NumWorkers = pbs_max_workers;

  job = createJob(cl);

  i = 1;
  for func = functions
    fprintf('Setting up job ID %d / %d...\n', i, length(functions));
    tasks(i) = createTask(job, @exp_saACMES_task, 0, {cwd, func, dimensions});
    i = i + 1;
  end

  tasks

  submit(job)
end
