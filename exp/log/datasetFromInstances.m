function dataset = datasetFromInstances(opts, nSnapshots, fun, dim, id)
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
  opts.maxEval = defopts(opts, 'maxEval', 250);
  
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

  % TODO:
  % cycle through instance saved in exp_settings.instances

  % BBOB fitness initialization
  fgeneric('initialize', exp_settings.bbob_function, exp_settings.instances(1), ['/tmp/bbob_output/']);

  % find maximal evaluation generation
  maxGener = find(cmaes_out{1}{1}.origEvaled, opts.maxEval*dim, 'first');
  nGenerations = cmaes_out{1}{1}.generations(maxGener(end));
  gens = floor(linspace(0, nGenerations-1, nSnapshots+1));
  gens(1) = [];
  cmo = cmaes_out{1}{1};

  % TODO:
  % make all these train/test sets 2-dimensional since we need all 15 instances
  % and move before line# 46
  trainSetX = cell(nSnapshots, 1);
  trainSetY = cell(nSnapshots, 1);
  testSetX = cell(nSnapshots, 1);
  testSetY = cell(nSnapshots, 1);
  
  % Dataset generation

  for i = 1:nSnapshots
    g = gens(i);

    lambda = sum(cmo.generations == g);

      % the model is not saved, so create a fresh new one
      xmean = cmo.means(:,g);
      sigma = cmo.sigmas(g);
      BD    = cmo.BDs{g};
      dim   = exp_settings.dim;
      m = ModelFactory.createModel(SF.surrogateParams.modelType, SF.surrogateParams.modelOpts, xmean');

      % create the Archive of original points
      minTrainSize = m.getNTrainData();
      archive = Archive(dim);
      orig_id = logical(cmo.origEvaled(1:(cmo.generationStarts(g)-1)));
      X_train = cmo.arxvalids(:,orig_id)';
      y_train = cmo.fvalues(orig_id)';
      archive.save(X_train, y_train, 5);
      archive.gens = cmo.generations(orig_id);

      % find the model's trainset -- points near the current xmean
      nArchivePoints = myeval(SF.surrogateParams.evoControlTrainNArchivePoints);
      [X_train, y_train, nData] = archive.getDataNearPoint(nArchivePoints, ...
        xmean', SF.surrogateParams.evoControlTrainRange, ...
        sigma, BD);
    end

    fprintf('Train set size (gen.# %3d): %3d\n', g, size(X_train,1));

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

      
    % set seed due to reproducibility of default dataset
    rng(opts.maxEval)
    
    % Generate fresh CMA-ES' \lambda offsprings
    [arx, arxvalid, arz] = sampleCmaesNoFitness(sigma, lambda, cmaesState, sampleOpts);

    % Save everything needed
    trainSetX{i} = X_train;
    trainSetY{i} = y_train;
    testSetX{i}  = arxvalid';
    testSetY{i}  = fgeneric(testSetX{i}')';
    means{i}     = xmean';
    sigmas{i}    = sigma;
    BDs{i}       = BD;
    cmaesStates{i} = cmaesState;
  end

  % Finalize

  dataset = struct();
  dataset.trainSetX = trainSetX;
  dataset.trainSetY = trainSetY;
  dataset.testSetX  = testSetX;
  dataset.testSetY  = testSetY;
  dataset.means = means;
  dataset.sigmas = sigmas;
  dataset.BDs = BDs;
  dataset.cmaesStates = cmaesStates;
  dataset.generations = gens;
  dataset.function  = fun;
  dataset.dim       = m.dim;
  dataset.id        = id;
  dataset.sampleOpts = sampleOpts;

  fgeneric('finalize');
end
