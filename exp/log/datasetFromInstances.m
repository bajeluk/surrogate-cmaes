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

  % sampleMethod:
  %   equidistant: deterministic with equal spacing
  %   uniform_wor: random sample without replacement, uniform
  %   weights
  %   uniform: random sample with replacement, uniform weights
  %   geometric: random sample with replacement, geometrical weights
  %   geometric: random sample with replacement, geometrical weights
  %   geometric_wor: random sample without replacement, geometrical weights
  opts.sampleMethod = defopts(opts, 'sampleMethod', 'equidistant');

  if ~exist('isForModelPool', 'var')
    isForModelPool = false;
    nPreviousGenerations = 0;
  end
  if (~exist('loadModels', 'var') || ~islogical(loadModels))
    loadModels = false;
  end
  % prepare output variable
  nInstances = length(inst);
  nIds = length(id);
  dataset = cell(nInstances, nIds);

  for ii_id = 1:nIds
    id_no = id(ii_id);

    % load data from S-CMA-ES log files
    scmaesOutFile = sprintf('%s/%s_results_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id_no);
    try
      SF = load(scmaesOutFile, 'cmaes_out', 'exp_settings', 'surrogateParams');
      cmaes_out = SF.cmaes_out;
      exp_settings = SF.exp_settings;
      exp_instances = SF.exp_settings.instances;
      if (loadModels)
        modellogOutFile = sprintf('%s/bbob_output/%s_modellog_%d_%dD_%d.mat', opts.exppath, exp_id, fun, dim, id_no);
        MF = load(modellogOutFile, 'models', 'models2');
      end
    catch
      warning('Model file or scmaes output file is missing in f%d %dD (id %d).', fun, dim, id_no)
      continue;
    end

    % cycle through instances
    for i_inst = 1:nInstances
      instanceNo = inst(i_inst);
      expInstanceId = find(exp_instances == instanceNo, 1);

      % look into the cmaes_out whether there this instance really is
      if (isempty(expInstanceId) || i_inst > length(SF.cmaes_out) || isempty(SF.cmaes_out{expInstanceId}{1}))
        warning('Instance %d is missing in the ''%s'' file.', instanceNo, scmaesOutFile);
        dataset{i_inst, ii_id} = [];
        % and skip this instance if not
        continue;
      end

      % BBOB fitness initialization
      user = getenv('USER');
      BBOB_PATH = ['/tmp/' user '/bbob_output'];
      [~, ~] = mkdir(BBOB_PATH);
      fgeneric('initialize', exp_settings.bbob_function, instanceNo, BBOB_PATH);

      % identify snapshot generations
      cmo = cmaes_out{expInstanceId}{1};
      % identify the first maxEval*dim orig-evaluated points
      origEvaledIdMax = find(cmo.origEvaled, opts.maxEval*dim, 'first');
      lastOrigEvaledIdMax = origEvaledIdMax(end);

      firstGeneration = cmo.generations(find(~cmo.origEvaled, 1));
      lastGeneration = cmo.generations(end);
      if (loadModels)
        lastGeneration = min(lastGeneration, length(MF.models));
      end
      nGenerations = lastGeneration-firstGeneration+1;

      if nSnapshots < 1
        % nSnapshots relative to all admissible generations
        nSnapshots = ceil(nSnapshots * nGenerations);
      end

      if (opts.uniqueGenerations || contains(opts.sampleMethod, 'wor'))
        % do not save the same generation multiple times, take only available
        % number of snapshots
        nSnapshots = min(nSnapshots, nGenerations);
        % gens = unique(gens);
      end

      % third: sample nSnapshots generations
      switch opts.sampleMethod
        case 'uniform_wor'
          gens = sort(randsample(firstGeneration:lastGeneration, nSnapshots, false));
        case 'uniform'
          gens = sort(randsample(firstGeneration:lastGeneration, nSnapshots, true));
        case 'geometric'
          % setting the parameter p of geometric distribution such that
          % nGenerations is the geometric distribution's quantile alpha
          alpha = 0.9;
          p = 1-nthroot(1-alpha, nGenerations+1);
          % setting weights by probability mass function of the geometric
          % distribution
          w = geopdf(0:nGenerations-1, p);
          gens = sort(randsample(firstGeneration:lastGeneration, nSnapshots, true, w));
        case 'geometric_wor'
          alpha = 0.9;
          p = 1-nthroot(1-alpha, nGenerations+1);
          w = geopdf(0:nGenerations-1, p);

          % drawing without replacement
          gens = [];
          draws = 0;
          while length(gens) < nSnapshots
            draws = draws + 1;
            draw = randsample(firstGeneration:lastGeneration, 1, true, w);

            i_inst, ii_id = find(gens >= draw, 1);
            if i_inst, ii_id && gens(i_inst, ii_id) > draw
              % insert sort with uniqueness
              gens = [gens(1:i_inst, ii_id-1) draw gens(i_inst, ii_id:end)];
            end

            if draws >= 1000 * nSnapshots
              warning('Could not draw a sample of required size without replacement, gave up.');
              break;
            end
          end
        case 'equidistant'
          gens = floor(linspace(firstGeneration, lastGeneration, nSnapshots));
        otherwise % case 'equidistant'
          warning('sampleMethod ''%s not'' recognized, fall back to ''equidistant''', ...
            opts.sampleMethod);
          gens = floor(linspace(firstGeneration, lastGeneration, nSnapshots));
      end

      % prepare output fields
      dataset{i_inst, ii_id} = struct();
      dataset{i_inst, ii_id}.archive  = [];
      if (isForModelPool)
        dataset{i_inst, ii_id}.testSetX = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.testSetY = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.means    = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.sigmas   = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.BDs      = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.diagDs   = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.diagCs   = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.pcs      = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.pss      = cell(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.iruns    = zeros(nSnapshots, nPreviousGenerations+1);
        dataset{i_inst, ii_id}.cmaesStates = cell(nSnapshots, nPreviousGenerations+1);
      else
        dataset{i_inst, ii_id}.testSetX = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.testSetY = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.means    = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.sigmas   = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.BDs      = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.diagDs   = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.diagDs   = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.pcs      = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.pss      = cell(1, nSnapshots);
        dataset{i_inst, ii_id}.iruns    = zeros(1, nSnapshots);
        dataset{i_inst, ii_id}.cmaesStates = cell(1, nSnapshots);
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
          diagD = cmo.diagDs{g};
          diagC = cmo.diagCs{g};
          pc    = cmo.pcs{g};
          ps    = cmo.pss{g};
          irun  = cmo.iruns(g);
          dim   = exp_settings.dim;

          cmaesState = struct( ...
            'xmean', xmean, ...
            'sigma', sigma, ...
            'lambda', lambda, ...
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
            dataset{i_inst, ii_id}.testSetX{sni, genShift+1}    = arxvalid';
            dataset{i_inst, ii_id}.testSetY{sni, genShift+1}    = fgeneric(arxvalid)';
            dataset{i_inst, ii_id}.means{sni, genShift+1}       = xmean';
            dataset{i_inst, ii_id}.sigmas{sni, genShift+1}      = sigma;
            dataset{i_inst, ii_id}.BDs{sni, genShift+1}         = BD;
            dataset{i_inst, ii_id}.cmaesStates{sni, genShift+1} = cmaesState;
            dataset{i_inst, ii_id}.sampleOpts{sni, genShift+1}  = sampleOpts;
            dataset{i_inst, ii_id}.diagDs{sni, genShift+1}      = diagD;
            dataset{i_inst, ii_id}.diagCs{sni, genShift+1}      = diagC;
            dataset{i_inst, ii_id}.pcs{sni, genShift+1}         = pc;
            dataset{i_inst, ii_id}.pss{sni, genShift+1}         = ps;
            dataset{i_inst, ii_id}.iruns(sni, genShift+1)       = irun;
          else
            dataset{i_inst, ii_id}.testSetX{sni}  = arxvalid';
            dataset{i_inst, ii_id}.testSetY{sni}  = fgeneric(arxvalid)';
            dataset{i_inst, ii_id}.means{sni}     = xmean';
            dataset{i_inst, ii_id}.sigmas{sni}    = sigma;
            dataset{i_inst, ii_id}.BDs{sni}       = BD;
            dataset{i_inst, ii_id}.cmaesStates{sni} = cmaesState;
            dataset{i_inst, ii_id}.sampleOpts{sni}  = sampleOpts;
            dataset{i_inst, ii_id}.diagDs{sni}      = diagD;
            dataset{i_inst, ii_id}.diagCs{sni}      = diagC;
            dataset{i_inst, ii_id}.pcs{sni}         = pc;
            dataset{i_inst, ii_id}.pss{sni}         = ps;
            dataset{i_inst, ii_id}.iruns(sni)       = irun;
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
                dataset{i_inst, ii_id}.models{sni}  = MF.models{thisG};
              else
                dataset{i_inst, ii_id}.models{sni}  = [];
              end
              % Try to save also the second model, if it is from that generation
              if (length(MF.models2) >= thisG && ~isempty(MF.models2{thisG}) ...
                  && MF.models2{thisG}.isTrained() && MF.models2{thisG}.trainGeneration == g)
                dataset{i_inst, ii_id}.models2{sni} = MF.models2{thisG};
              else
                dataset{i_inst, ii_id}.models2{sni}  = [];
              end
            else
              dataset{i_inst, ii_id}.models{sni}  = [];
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

      dataset{i_inst, ii_id}.archive     = archive;
      dataset{i_inst, ii_id}.generations = gens;
      dataset{i_inst, ii_id}.function  = fun;
      dataset{i_inst, ii_id}.dim       = dim;
      dataset{i_inst, ii_id}.id        = id_no;
      dataset{i_inst, ii_id}.instance  = instanceNo;
%       dataset{i_inst, ii_id}.maxEval   = opts.maxEval;

      fgeneric('finalize');
    end  % instances loop
  end % ids loop
end
