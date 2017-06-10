function dataset = datasetFromInstances(opts, nSnapshots, fun, dim, inst, id, isForModelPool, nPreviousGenerations, loadModels)
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

  % start processing result on generation where this number of FE/D is
  % achieved
  opts.startEval = defopts(opts, 'startEval', 1);
  % process the results only until this number of evalutions per dimension
  opts.maxEval = defopts(opts, 'maxEval', 250);
  % make multiple testcases per generation even if there is more snapshots 
  % than generations to test
  opts.uniqueGenerations = defopts(opts, 'uniqueGenerations', false);

  if ~exist('isForModelPool', 'var')
    isForModelPool = false;
    nPreviousGenerations = 0;
  end
  if (~exist('loadModels', 'var') || ~islogical(loadModels))
    loadModels = false;
  end
  % prepare output variable
  nInstances = length(inst);
  dataset = cell(1, nInstances);

  % load data from S-CMA-ES log files
  scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id);
  try
    SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
    cmaes_out = SF.cmaes_out;
    exp_settings = SF.exp_settings;
    exp_instances = SF.exp_settings.instances;
    if (loadModels)
      modellogOutFile = sprintf('%s/bbob_output/%s_modellog_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id);
      MF = load(modellogOutFile, 'models', 'models2');
    end
  catch
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

    % BBOB fitness initialization
    user = getenv('USER');
    BBOB_PATH = ['/tmp/' user '/bbob_output'];
    [~, ~] = mkdir(BBOB_PATH);
    fgeneric('initialize', exp_settings.bbob_function, instanceNo, BBOB_PATH);

    % identify snapshot generations
    % TODO: make exponential gaps between snapshot generations
    %       to have higher density at start of optim. run
    %
    cmo = cmaes_out{expInstanceId}{1};
    % first, identify the first maxEval*dim orig-evaluated points
    origEvaledIdStart = find(cmo.origEvaled, opts.startEval*dim, 'first');
    lastOrigEvaledIdStart = origEvaledIdStart(end);
    origEvaledIdMax = find(cmo.origEvaled, opts.maxEval*dim, 'first');
    lastOrigEvaledIdMax = origEvaledIdMax(end);
    % second, identify its generation
    firstGeneration = cmo.generations(lastOrigEvaledIdStart);
    lastGeneration = cmo.generations(lastOrigEvaledIdMax);
    if (loadModels)
      lastGeneration = min(lastGeneration, length(MF.models));
    end
    if (opts.uniqueGenerations)
      % do not save the same generation multiple times, take only available
      % number of snapshots
      nSnapshots = min(nSnapshots, lastGeneration - firstGeneration + 1);
      % gens = unique(gens);
    end
    % third place the snapshot generations equdistant in generations...
    gens = floor(linspace(firstGeneration, lastGeneration, nSnapshots));

    % prepare output fields
    dataset{i_inst} = struct();
    dataset{i_inst}.archive  = [];
    if (isForModelPool)
      dataset{i_inst}.testSetX = cell(nSnapshots, nPreviousGenerations+1);
      dataset{i_inst}.testSetY = cell(nSnapshots, nPreviousGenerations+1);
      dataset{i_inst}.means    = cell(nSnapshots, nPreviousGenerations+1);
      dataset{i_inst}.sigmas   = cell(nSnapshots, nPreviousGenerations+1);
      dataset{i_inst}.BDs      = cell(nSnapshots, nPreviousGenerations+1);
      dataset{i_inst}.cmaesStates = cell(nSnapshots, nPreviousGenerations+1);
    else
      dataset{i_inst}.testSetX = cell(1, nSnapshots);
      dataset{i_inst}.testSetY = cell(1, nSnapshots);
      dataset{i_inst}.means    = cell(1, nSnapshots);
      dataset{i_inst}.sigmas   = cell(1, nSnapshots);
      dataset{i_inst}.BDs      = cell(1, nSnapshots);
      dataset{i_inst}.cmaesStates = cell(1, nSnapshots);
    end

    % Dataset generation

    for sni = 1:nSnapshots        % sni stands for SNapshot Index
      for genShift = nPreviousGenerations:-1:0
        g = gens(sni);
        g = g-genShift;
        if (g < 1)
          continue;
        end

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

        if (isfield(cmo, 'arxvalids') && ~isempty(cmo.arxvalids))
          arxvalid = cmo.arxvalids(:, cmo.generations == g);
        else
          % Generate fresh CMA-ES' \lambda offsprings
          [~, arxvalid, ~] = sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
        end

        % Save everything needed
        if (isForModelPool)
          dataset{i_inst}.testSetX{sni, genShift+1}    = arxvalid';
          dataset{i_inst}.testSetY{sni, genShift+1}    = fgeneric(arxvalid)';
          dataset{i_inst}.means{sni, genShift+1}       = xmean';
          dataset{i_inst}.sigmas{sni, genShift+1}      = sigma;
          dataset{i_inst}.BDs{sni, genShift+1}         = BD;
          dataset{i_inst}.cmaesStates{sni, genShift+1} = cmaesState;
          dataset{i_inst}.sampleOpts{sni, genShift+1}  = sampleOpts;
        else
          dataset{i_inst}.testSetX{sni}  = arxvalid';
          dataset{i_inst}.testSetY{sni}  = fgeneric(arxvalid)';
          dataset{i_inst}.means{sni}     = xmean';
          dataset{i_inst}.sigmas{sni}    = sigma;
          dataset{i_inst}.BDs{sni}       = BD;
          dataset{i_inst}.cmaesStates{sni} = cmaesState;
          dataset{i_inst}.sampleOpts{sni} = sampleOpts;
        end

        % save models if required
        if (loadModels)
          % save only models from within one generation from current generation 'g'
          if (length(MF.models) >= g && ~isempty(MF.models{g}) ...
              && MF.models{g}.isTrained() && abs(MF.models{g}.trainGeneration - g) <= 1)
            % Sometimes it happens that we have saved a model one generation younger...
            if (MF.models{g}.trainGeneration - g == 1)
              thisG = g-1;
            else
              thisG = g;
            end
            if (MF.models{thisG}.isTrained())
              dataset{i_inst}.models{sni}  = MF.models{thisG};
            else
              dataset{i_inst}.models{sni}  = [];
            end
            % Try to save also the second model, if it is from that generation
            if (length(MF.models2) >= thisG && ~isempty(MF.models2{thisG}) ...
                && MF.models2{thisG}.isTrained() && MF.models2{thisG}.trainGeneration == g)
              dataset{i_inst}.models2{sni} = MF.models2{thisG};
            else
              dataset{i_inst}.models2{sni}  = [];
            end
          else
            dataset{i_inst}.models{sni}  = [];
          end
        end

      end % generationShift loop
    end  % snapshots loop

    % create the Archive of with the original-evaluated points
    archive = Archive(dim);
    orig_id = logical(cmo.origEvaled(1:lastOrigEvaledIdMax));
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
