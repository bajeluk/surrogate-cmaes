function dataset = datasetFromInstances(opts, nSnapshots, fun, dim, inst, id)
%DATASETFROMINSTANCE - generates datasets for specified dim and #fun for offline model tunning

% Generates datasets for offline model tunning from the DTS-CMA-ES
% results files with 'progressLog' switched on (*_results_*.mat)

  % Initialization

  % DEBUG:
  % exp_id = 'exp_doubleEC_21_log';
  % nSnapshots = 10;
  % fun = 1;
  % dim = 2;
  % id = 1;

  if nargin < 1
    help datasetFromInstance
    return
  end
  exp_id = opts.inputExp_id;

  opts.maxEval = defopts(opts, 'maxEval', 250);

  % prepare output variable
  nInstances = length(inst);
  dataset = cell(1, nInstances);

  % load data from S-CMA-ES log files
  scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id);
  if exist(scmaesOutFile, 'file')
    SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
    cmaes_out = SF.cmaes_out;
    exp_settings = SF.exp_settings;
    exp_instances = SF.exp_settings.instances;
  else
    warning('Model file or scmaes output file is missing in f%d %dD (id %d).', fun, dim, id)
    return
  end

  % cycle through instances
  for i_inst = 1:nInstances
    instanceNo = inst(i_inst);
    expInstanceId = find(exp_instances == instanceNo, 1);

    % look into the cmaes_out whether there this instance really is
    if (isempty(expInstanceId) || i_inst > length(SF.cmaes_out) || isempty(SF.cmaes_out{expInstanceId}{1}))
      warning('Instance %d is missing in the ''%s'' file.', instanceNo, scmaesOutFile);
      dataset{i_inst} = [];
      % and skip this instance if not
      continue;
    end

    % prepare output fields
    dataset{i_inst} = struct();
    dataset{i_inst}.archive  = [];
    dataset{i_inst}.testSetX = cell(1, nSnapshots);
    dataset{i_inst}.testSetY = cell(1, nSnapshots);
    dataset{i_inst}.means    = cell(1, nSnapshots);
    dataset{i_inst}.sigmas   = cell(1, nSnapshots);
    dataset{i_inst}.BDs      = cell(1, nSnapshots);
    dataset{i_inst}.cmaesStates = cell(1, nSnapshots);
    dataset{i_inst}.generations = [];

    % BBOB fitness initialization
    fgeneric('initialize', exp_settings.bbob_function, instanceNo, '/tmp/bbob_output/');

    % identify snapshot generations
    % TODO: make exponential gaps between snapshot generations
    %       to have higher density at start of optim. run
    %
    cmo = cmaes_out{expInstanceId}{1};
    % first, identify the first maxEval*dim orig-evaluated points
    origEvaledId = find(cmo.origEvaled, opts.maxEval*dim, 'first');
    lastOrigEvaledId = origEvaledId(end);
    % second, identify its generation
    lastGeneration = cmo.generations(lastOrigEvaledId);
    % third place the snapshot generations equdistant in generations...
    gens = floor(linspace(2, lastGeneration, nSnapshots+1));
    % ...but do not start at the beginning
    gens(1) = [];

    % Dataset generation

    for sni = 1:nSnapshots        % sni stands for SNapshot Index
      g = gens(sni);

      lambda = sum(cmo.generations == g);

      xmean = cmo.means(:,g);
      sigma = cmo.sigmas(g);
      BD    = cmo.BDs{g};
      dim   = exp_settings.dim;

      cmaesState = struct( ...
        'xmean', xmean, ...
        'sigma', sigma, ...
        'lambda', lambda, ...
        'BD', BD, ...
        ... % 'diagD', diagD, ...
        'diagD', [], ...
        ... % 'diagC', diagC, ...
        'dim', dim, ...
        'mu', floor(lambda/2), ...
        'countiter', g);

      sampleOpts = struct( ...
        'noiseReevals', 0, ...
        'isBoundActive', true, ...
        'lbounds', -5 * ones(dim, 1), ...
        'ubounds',  5 * ones(dim, 1), ...
        'counteval', cmo.generationStarts(g), ...
        'flgEvalParallel', false, ...
        'flgDiagonalOnly', false, ...
        'noiseHandling', false, ...
        'xintobounds', @xintobounds, ...
        'origPopSize', lambda);

      % Generate fresh CMA-ES' \lambda offsprings
      [~, arxvalid, ~] = sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);

      % Save everything needed
      dataset{i_inst}.testSetX{sni}  = arxvalid';
      dataset{i_inst}.testSetY{sni}  = fgeneric(arxvalid)';
      dataset{i_inst}.means{sni}     = xmean';
      dataset{i_inst}.sigmas{sni}    = sigma;
      dataset{i_inst}.BDs{sni}       = BD;
      dataset{i_inst}.cmaesStates{sni} = cmaesState;
      dataset{i_inst}.sampleOpts{sni} = sampleOpts;

    end  % snapshots loop

    % create the Archive of with the original-evaluated points
    archive = Archive(dim);
    orig_id = logical(cmo.origEvaled(1:lastOrigEvaledId));
    X_orig = cmo.arxvalids(:,orig_id)';
    y_orig = cmo.fvalues(orig_id)';
    archive.save(X_orig, y_orig, cmo.generations(orig_id));

    dataset{i_inst}.archive     = archive;
    dataset{i_inst}.generations = gens;
    dataset{i_inst}.function  = fun;
    dataset{i_inst}.dim       = dim;
    dataset{i_inst}.id        = id;
    dataset{i_inst}.instance  = instanceNo;
    dataset{i_inst}.maxEval   = opts.maxEval;

    fgeneric('finalize');
  end  % instances loop
end
