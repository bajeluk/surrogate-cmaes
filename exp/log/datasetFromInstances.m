function dataset = datasetFromInstances(exp_id, fun, dim, inst, id, varargin)
%DATASETFROMINSTANCES - generates datasets for specified dim and #fun for 
% offline model tunning
% 
% dataset = datasetFromInstances(exp_id, fun, dim, inst, id, opts)
% Generates datasets for offline model tunning from the DTS-CMA-ES
% results files with 'progressLog' switched on (*_results_*.mat).
%
% Input:
%   exp_id - experiment id where the CMA-ES logs are saved | string |
%            default: 'exp_doubleEC_21_log'
%   fun    - BBOB function number to load | positive integer scalar |
%            default: 1
%   dim    - dimension to load (2, 3, 5, 10, ...) | positive integer scalar
%            | default: 2
%   inst   - instance to load (1:5, ...) | positive integer scalar |
%            default: 1
%   id     - job id to load (1, 2, ..., # of possible combinations) |
%            positive integer scalar | default: 1

%   opts   - pairs of property (string) and value or struct with properties
%            as fields:
%     'expPath'              - path to experiment data folder | string |
%                              default: exp/experiments/exp_id
%     'isForFeatures'        - dataset for features testing | boolean |
%                              default: false
%     'isForModelPool'       - dataset for ModelPool testing | boolean |
%                              default: false
%     'loadModels'           - load saved models for ModelPool | boolean |
%                              default: false
%     'maxEval'              - process the results only until this number
%                              of evalutions per dimension | positive
%                              integer scalar | default: 250
%     'nPreviousGenerations' - number of previous generations to load for
%                              ModelPool testing | integer scalar |
%                              default: 0
%     'nSampleArchives'      - number of generated archives per one CMA-ES
%                              run for feature testing | positive integer
%                              scalar | default: 10
%     'nSmoothGenerations'   - number of extra generations for distribution
%                              smoothing in feature testing (window radius)
%                              | non-negative integer scalar | default: 0
%                            - example: nSmoothGenerations = 2 => window
%                              size = 2 (previous gens) + 1 (actual gen.) +
%                              2 (following gens) = 5 gens
%     'nSnapshotsPerRun'     - number of generated datasets per one CMA-ES
%                              run | positive integer scalar | default: 10
%     'sampleMethod'         - generation sample method | string (see
%                              below) | default: 'equidistant'
%        'equidistant'   - deterministic with equal spacing
%        'geometric'     - random sample with replacement, geometrical
%                          weights
%        'geometric_wor' - random sample without replacement, geometrical
%                          weights
%        'uniform'       - random sample with replacement, uniform weights
%        'uniform_wor'   - random sample without replacement, uniform
%                          weights
%     'smoothDistribution'   - floating window style of distribution
%                              weighting for feature testing | string:
%                              'average', 'exponential', 'gaussian',
%                              'identical', 'linear' | default: 'gaussian'
%     'startEval'            - start processing result on generation where
%                              this number of FE/D is achieved | positive
%                              integer scalar | default: 1
%     'uniqueGenerations'    - make multiple testcases per generation even
%                              if there is more snapshots than generations
%                              to test | boolean | default: false
%
% Output:
%   dataset - #inst x #id cell-array of structures with the DTS-CMA-ES
%             state variables as fields
%
% See Also:
%   modelTestSets

  % Initialization
  if nargout > 0
    dataset = {};
  end

  if nargin < 5
    if nargin < 4
      if nargin < 3
        if nargin < 2
          if nargin < 1
            help datasetFromInstances
            return
          end
          fun = 1;
        end
        dim = 2;
      end
      inst = 1;
    end
    id = 1;
  end

  % parse settings
  opts = settings2struct(varargin{:});

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
  %   geometric_wor: random sample without replacement, geometrical weights
  opts.sampleMethod = defopts(opts, 'sampleMethod', 'equidistant');
  % path to experiment folder
  opts.exppath = defoptsi(opts, 'exppath', fullfile('exp', 'experiments', exp_id));

  % model pool settings
  isForModelPool = defopts(opts, 'isForModelPool', false);
  nPreviousGenerations = defopts(opts, 'nPreviousGenerations', 0);
  loadModels = defopts(opts, 'loadModels', false);

  % feature testing settings
  isForFeatures = defopts(opts, 'isForFeatures', false);
  nSmoothGenerations = defopts(opts, 'nSmoothGenerations', 0);
  smoothDistribution = defopts(opts, 'smoothDistribution', 'gaussian');
  nSampleArchives = defopts(opts, 'nSampleArchives', 10);

  % number of generation snapshots
  nSnapshots = defopts(opts, 'nSnapshotsPerRun', 10);

  % prepare output variable
  nInstances = length(inst);
  nIds = length(id);
  dataset = cell(nInstances, nIds);

  % id (model settings) cycle
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
      % identify all the original-evaluated points
      orig_id = find(cmo.origEvaled);
      % identify unique original-evaluated points
      [~, uniqueOrigId] = unique(cmo.arxvalids(:, cmo.origEvaled)', 'rows');
      % select ids of unique points preserving the original ordering
      orig_id = orig_id(sort(uniqueOrigId));
      % select maxEval*dim points at maximum
      orig_id = orig_id(1:min(numel(orig_id), opts.maxEval*dim));

      firstModelGeneration = cmo.generations(find(~cmo.origEvaled, 1));
      lastGeneration = cmo.generations(end);
      if (loadModels)
        lastGeneration = min(lastGeneration, length(MF.models));
      end
      nGenerations = lastGeneration - firstModelGeneration + 1;

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
          gens = sort(randsample(firstModelGeneration:lastGeneration, nSnapshots, false));
        case 'uniform'
          gens = sort(randsample(firstModelGeneration:lastGeneration, nSnapshots, true));
        case 'geometric'
          % setting the parameter p of geometric distribution such that
          % nGenerations is the geometric distribution's quantile alpha
          alpha = 0.9;
          p = 1-nthroot(1-alpha, nGenerations+1);
          % setting weights by probability mass function of the geometric
          % distribution
          w = geopdf(0:nGenerations-1, p);
          gens = sort(randsample(firstModelGeneration:lastGeneration, nSnapshots, true, w));
        case 'geometric_wor'
          % TODO: make this part of code more readable
          alpha = 0.9;
          p = 1-nthroot(1-alpha, nGenerations+1);
          w = geopdf(0:nGenerations-1, p);

          % drawing without replacement
          gens = [];
          draws = 0;
          while length(gens) < nSnapshots
            draws = draws + 1;
            draw = randsample(firstModelGeneration:lastGeneration, 1, true, w);

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
          gens = floor(linspace(firstModelGeneration, lastGeneration, nSnapshots));
        otherwise % case 'equidistant'
          warning('sampleMethod ''%s not'' recognized, fall back to ''equidistant''', ...
            opts.sampleMethod);
          gens = floor(linspace(firstModelGeneration, lastGeneration, nSnapshots));
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
      elseif (isForFeatures)
        dataset{i_inst, ii_id}.testSetX = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.testSetY = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.means    = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.sigmas   = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.BDs      = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.diagDs   = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.diagCs   = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.pcs      = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.pss      = cell(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.iruns    = zeros(nSnapshots, nSampleArchives);
        dataset{i_inst, ii_id}.cmaesStates = cell(nSnapshots, nSampleArchives);
      else
        dataset{i_inst, ii_id}.testSetX = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.testSetY = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.means    = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.sigmas   = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.BDs      = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.diagDs   = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.diagDs   = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.pcs      = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.pss      = cell(nSnapshots, 1);
        dataset{i_inst, ii_id}.iruns    = zeros(nSnapshots, 1);
        dataset{i_inst, ii_id}.cmaesStates = cell(nSnapshots, 1);
      end

      % Dataset generation

      for sni = 1:nSnapshots        % sni stands for SNapshot Index
        for genShift = nPreviousGenerations:-1:0
          % original generation number
          g = gens(sni);
          g = g-genShift;
          if (g < 1)
            continue;
          end

          % CMA-ES state variables
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
            'BD', BD, ...
            'diagD', diagD, ...
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
          elseif (isForFeatures)
            for sai = 1:nSampleArchives % sai = Sample Archive Index
              % Generate fresh CMA-ES' \lambda offsprings for each sample
              [~, arxvalid, ~] = sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);
              dataset{i_inst, ii_id}.testSetX{sni, sai}    = arxvalid';
              dataset{i_inst, ii_id}.testSetY{sni, sai}    = fgeneric(arxvalid)';
            end
            % means, sigmas, BDs, cmaesStates, sampleOpts, diagDs, and
            % diagCs will be generated later during Archives generation
            % process
            dataset{i_inst, ii_id}.pcs(sni, :)         = {pc};
            dataset{i_inst, ii_id}.pss(sni, :)         = {ps};
            dataset{i_inst, ii_id}.iruns(sni, :)       = irun;
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

      if (isForFeatures)
        % create Archives with points sampled from smoothed CMA-ES
        % distribution through a range of generations
        archive = cell(1, nSampleArchives);
        for sai = 1:nSampleArchives % sai = Sample Archive Index
          archive{sai} = Archive(dim);
        end
        generations = cmo.generations(orig_id);
        % save initial point to individual archives
        if ismember(0, generations)
          for sai = 1:nSampleArchives
            % save the first mean as starting point (as it is in CMA-ES)
            archive{sai}.save(cmo.means(:, 1)', fgeneric(cmo.means(:, 1)), 0);
          end
        end
        % archive generating through distribution smoothing
        for gi = 1:lastGeneration
          % sample for all archives at once
          [X_smooth, cmaesState, sampleOpts] = smoothedSampling(...
            cmo, gi, nSampleArchives, smoothDistribution, nSmoothGenerations);

          % save points in generations with at least one original evaluated
          % point
          if ismember(gi, generations)
            % evaluate points for all archives
            y_smooth = fgeneric(X_smooth)';
            % point identifier in Xsmooth
            pId = 0;
            % save points to individual archives
            for sai = 1:nSampleArchives
              % while there are points to save to the archive in current
              % generation
              % TODO: this cycle can be slow, on the other hand, ensures
              %       the same numbers of points in archives regardless
              %       duplicities
              while numel(archive{sai}.y) < sum(generations <= gi)
                % get remaining ids of generated points
                pId = pId(end) + (1:(sum(generations <= gi) - numel(archive{sai}.y)));
                % in case of reaching the last point, generate new points
                if any(pId > numel(y_smooth))
                  pId = 1:(sum(generations <= gi) - numel(archive{sai}.y));
                  % Generate fresh CMA-ES' offspring
                  [~, X_smooth, ~] = sampleCmaesNoFitness( ...
                    cmaesState.sigma, nSampleArchives*cmaesState.lambda, ...
                    cmaesState, sampleOpts);
                  y_smooth = fgeneric(X_smooth)';
                end
                % save original evaluated points
                archive{sai}.save(X_smooth(:, pId)', y_smooth(pId), gi);
              end
            end
          end

          % save cmaesState if gi == snapshot generation
          snapshotGen = (gi == gens);
          if any(snapshotGen)
            dataset{i_inst, ii_id}.means(snapshotGen, :)       = {cmaesState.xmean'};
            dataset{i_inst, ii_id}.sigmas(snapshotGen, :)      = {cmaesState.sigma};
            dataset{i_inst, ii_id}.BDs(snapshotGen, :)         = {cmaesState.BD};
            dataset{i_inst, ii_id}.cmaesStates(snapshotGen, :) = {cmaesState};
            dataset{i_inst, ii_id}.sampleOpts(snapshotGen, :)  = {sampleOpts};
            dataset{i_inst, ii_id}.diagDs(snapshotGen, :)      = {cmaesState.diagD};
            dataset{i_inst, ii_id}.diagCs(snapshotGen, :)      = {cmaesState.diagC};
          end
        end
      else
        % create the Archive with the original-evaluated points
        archive = Archive(dim);
        X_orig = cmo.arxvalids(:,orig_id)';
        if isfield(cmo, 'fvaluesOrig')
          y_orig = cmo.fvaluesOrig(orig_id)';
        else
          y_orig = fgeneric(X_orig);
        end
        archive.save(X_orig, y_orig, cmo.generations(orig_id));
      end

      % save archive
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

function [arxvalid, cmaesState, sampleOpts] = smoothedSampling(cmo, gen, nPopulations, weightStyle, windowGens)
% SMOOTHEDSAMPLING - generate archive points in the input space using
% smoothing of the original CMA-ES distributions
%
% [arxvalid, cmaesState, sampleOpts] = ...
%   smoothedSampling(cmo, gen, nPopulations, weightStyle, windowGens)
%
% Input:
%   cmo          - cmaes out structure | structure
%   gen          - actual generation | positive integer scalar
%   nPopulations - number of populations to generate | positive integer
%                  scalar
%   weightStyle  - floating window style of distribution weighting |
%                  'average', 'exponential', 'gaussian', 'identity',
%                  'linear'
%   windowGens   - number of generations added to window range |
%                  non-negative integer scalar
%
% Output:
%   arxvalid   - sampled populations from smoothed CMA-ES distribution |
%                double matrix
%   cmaesState - CMA-ES state variables | structure
%   sampleOpts - sampling options | structure

  [dim, nGens] = size(cmo.means);
  % get weights
  weight = getWeights(gen, nGens, weightStyle, windowGens);

  % mean
  xmean = weight*cmo.means';
  % sigma
  sigma = weight*cmo.sigmas';
  % sigma*BD = \sum_{i=1}^{g_max} weight_i*sigma_i*BD_i
  sBDcell = arrayfun(@(x) weight(x)*cmo.sigmas(x)*cmo.BDs{x}, ...
              find(weight > 0), 'UniformOutput', false);
  % sum weighted sBD matrices
  BD = sum( cat(3, sBDcell{:}), 3) / sigma;
  % diagonal
  diagC = diag(BD*BD');

  % the same with diagonal covariance
  sdiagDcell = arrayfun(@(x) weight(x)*cmo.sigmas(x)*cmo.diagDs{x}, ...
               find(weight > 0), 'UniformOutput', false);
  diagD = sum( cat(3, sdiagDcell{:}), 3) / sigma;

  % set the rest of state variables
  lambda = sum(cmo.generations == gen);

  cmaesState = struct( ...
    'xmean', xmean', ...
    'sigma', sigma, ...
    'lambda', lambda, ...
    'dim', numel(xmean), ...
    'mu', floor(lambda/2), ...
    'BD', BD, ...
    'diagD', diagD, ...
    'diagC', diagC, ...
    'countiter', gen ...
  );

  sampleOpts = struct( ...
    'noiseReevals', 0, ...
    'isBoundActive', true, ...
    'lbounds', -5 * ones(dim, 1), ...
    'ubounds',  5 * ones(dim, 1), ...
    'counteval', cmo.generationStarts(gen), ...
    'flgEvalParallel', false, ...
    'flgDiagonalOnly', false, ...
    'noiseHandling', false, ...
    'xintobounds', @xintobounds, ...
    'origPopSize', lambda ...
  );

  % Generate fresh CMA-ES' offsprings
  [~, arxvalid, ~] = sampleCmaesNoFitness(sigma, nPopulations*lambda, cmaesState, sampleOpts);
end

function weights = getWeights(gen, nGens, weightStyle, windowGens)
% GETWEIGHTS - Generate weights for CMA-ES generation smoothing according
% to weightStyle and windowGens.
%
% Input:
%   gen         - current generation
%   nGens       - overall number of generations in particular run
%   weightStyle - style of weighting | 'average', 'exponential',
%                 'gaussian', 'identity', 'linear'
%
%   windowGens  - number of generations added to window range:
%                   0 - only original generation
%                       [..., 0, g, 0, ...]
%                   1 - one generation range
%                       [..., 0, g-1, g, g+1, 0, ...]
%                   2 - two generation range
%                       [..., 0, g-2, g-1, g, g+1, g-2, 0, ...]
%
% Output:
%   weights - weight vector

  if nargin < 4
    windowGens = 2;
  end

  weights = zeros(1, nGens);
  switch weightStyle
    case 'average'
      weights(:) = 1;
    case 'exponential'
      windowBase = exp([1:windowGens+1, windowGens:-1:1]);
    case 'gaussian'
      % normal distribution
      pd = makedist('Normal', 'mu', 0, 'sigma', 1);
      windowBase = pdf(pd, -windowGens:windowGens);
    case 'identity'
      weights(gen) = 1;
    case 'linear'
      windowBase = [1:windowGens+1, windowGens:-1:1];
    otherwise
      warning('There is no weighting style ''%s''. Setting ''identity''', weightStyle)
      weights(gen) = 1;
  end

  % window settings to weights
  if any(strcmp(weightStyle, {'exponential', 'gaussian', 'linear'}))
    % weight ids
    lowId = max(1, gen-windowGens);
    highId = min(nGens, gen+windowGens);
    % window ids
    winLowId = max(1, windowGens-gen+2);
    winHighId = min(windowGens, nGens-gen) + windowGens + 1;
    % add appropriate part of window
    weights(lowId:highId) = windowBase(winLowId:winHighId);
  end

  % normalize to unit lenght
  weights = weights/sum(weights);

end
