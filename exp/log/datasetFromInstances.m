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

  dataset = {};

  if nargin < 1
    help datasetFromInstance
    return
  end
  exp_id = opts.inputExp_id;
  
  % load data from files
  scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id);
  if exist(scmaesOutFile, 'file')
    SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
    cmaes_out = SF.cmaes_out;
    exp_settings = SF.exp_settings;
  else
    warning('Model file or scmaes output file is missing in f%d %dD (id %d).', fun, dim, id)
    return
  end

  opts.maxEval = defopts(opts, 'maxEval', 250);
  availInstances = intersect(inst, exp_settings.instances);
  if (length(availInstances) < length(inst))
    warning(['There are instances missing in the ' scmaesOutFile ' file.']);
  end
  opts.nInstances = length(availInstances);

  % prepare output variables
  dataset = struct();
  dataset.archives = cell(opts.nInstances, 1);
  dataset.testSetX = cell(opts.nInstances, nSnapshots);
  dataset.testSetY = cell(opts.nInstances, nSnapshots);
  dataset.means    = cell(opts.nInstances, nSnapshots);
  dataset.sigmas   = cell(opts.nInstances, nSnapshots);
  dataset.BDs      = cell(opts.nInstances, nSnapshots);
  dataset.cmaesStates = cell(opts.nInstances, nSnapshots);
  dataset.generations = cell(opts.nInstances, 1);

  for i_inst = 1:opts.nInstances
    instanceNo = availInstances;

    % BBOB fitness initialization
    fgeneric('initialize', exp_settings.bbob_function, instanceNo, '/tmp/bbob_output/');

    % identify snapshot generations
    %
    % first, identify index the first orig-evaluated point just exceeding 'maxEval'
    maxEvalGeneration = find(cmaes_out{1}{1}.origEvaled, opts.maxEval*dim, 'first');
    % second, identify its generation
    lastGeneration = cmaes_out{1}{1}.generations(maxEvalGeneration(end)) - 1;
    % third place the spapshot generations equdistant in generations
    gens = floor(linspace(0, lastGeneration, nSnapshots+1));
    gens(1) = [];
    cmo = cmaes_out{1}{1};

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
      dataset.testSetX{i_inst, sni}  = arxvalid';
      dataset.testSetY{i_inst, sni}  = fgeneric(dataset.testSetX{i_inst, sni}')';
      dataset.means{i_inst, sni}     = xmean';
      dataset.sigmas{i_inst, sni}    = sigma;
      dataset.BDs{i_inst, sni}       = BD;
      dataset.cmaesStates{i_inst, sni} = cmaesState;

    end  % snapshots loop

    % create the Archive of with the original-evaluated points
    archive = Archive(dim);
    orig_id = logical(cmo.origEvaled(1:lastGeneration));
    X_orig = cmo.arxvalids(:,orig_id)';
    y_orig = cmo.fvalues(orig_id)';
    archive.save(X_orig, y_orig, cmo.generations(orig_id));

    dataset.archives{i_inst}    = archive;
    dataset.generations{i_inst} = gens;

    fgeneric('finalize');
  end  % instances loop

  dataset.function  = fun;
  dataset.dim       = dim;
  dataset.id        = id;
  dataset.sampleOpts = sampleOpts;
  dataset.maxEval   = opts.maxEval;
end
