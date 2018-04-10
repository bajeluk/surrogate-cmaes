function testMetaLearn(modelOptions, modelOptionsInd, opts, funcs, dims, ...
  insts, Ns, designTypes)
%TESTMETALEARN Test models on metalearning data sets.

  if ~iscell(Ns), Ns = {Ns}; end
  if ~iscell(designTypes), designTypes = {designTypes}; end

  assert(isnumeric(funcs), '''funcs'' has to be numeric');
  assert(isnumeric(dims), '''dims'' has to be numeric');
  assert(isnumeric(insts), '''insts'' has to be numeric');

  % Default options
  opts.modelTypes = defopts(opts, 'modelTypes', {'rf', 'gp'});
  opts.cv_type = defopts(opts, 'cv_type', 'KFold');
  opts.cv_param = defopts(opts, 'cv_param', 10);
  opts.cv_ind = defopts(opts, 'cv_ind', 1:opts.cv_param);
  opts.inst_subsamp = defopts(opts, 'inst_subsamp', false);
  opts.use_rng_seed = defopts(opts, 'use_rng_seed', true);
  opts.rng_seed = defopts(opts, 'rng_seed', 42);
  opts.use_parpool = defopts(opts, 'use_parpool', true);
  opts.parpool_size = defopts(opts, 'parpool_size', 10);
  opts.dataset_path = defopts(opts, 'dataset_path', 'data_metalearning');
  opts.fname_template = defopts(opts, 'fname_template', strjoin({'data_metalearn_', ...
    '%dD_', ...
    'f%d_', ...
    'inst%d_', ...
    'N%d_', ...
    'design-%s.mat'}, '' ...
  ));
  opts.res_path = defopts(opts, 'res_path', 'results');
  opts.res_fname_template = defopts(opts, 'res_fname_template', strjoin({'res_', ...
    '%dD_', ...
    'f%d_', ...
    'inst%d_', ...
    'N%d_', ...
    'design-%s_', ...
    'model-%s_', ...
    'opts%d.mat'}, '' ...
  ));
  opts.rewrite_results = defopts(opts, 'rewrite_results', false);
  opts.save_retries = defopts(opts, 'save_retries', 10);

  % parpool init
%  if opts.use_parpool
%    w = MetaParPool('open');
%    fprintf('Started metacentrum pool with %d workers.\n', w);
%  end

  % dimension loop
  for dim = dims
    % function loop
    for func = funcs

      if opts.inst_subsamp
        % 0 in result name marks that instances are sampled in each CV fold
        insts_all = [0];
      else
        insts_all = insts;
      end

      % instance loop
      for inst = insts_all

        % dataset sizes loop
        for N_cell = Ns
          N = myeval(N_cell{:});

          % CV partitioning
          cv = make_cvpartition(opts, N);

          % design types loop
          for design_cell = designTypes
            designType = design_cell{:};

            % load the dataset
            if opts.inst_subsamp
              % instance subsampling
              inst_samp = randsample(insts, cv.NumTestSets, true); % with replacement
            else
              inst_samp = [inst];
            end
            fprintf('Inst sample: %s\n', num2str(inst_samp));

            data = cell(1, max(insts));
            for instance = inst_samp
              fname = sprintf(opts.fname_template, dim, func, instance, N, designType);
              fname = fullfile(opts.dataset_path, sprintf('%dD', dim), fname);
              d = load(fname);

              % data sanity checks
              assert(dim == size(d.X, 1), 'Unexpected dimensionality.');
              assert(N == size(d.X, 2), 'Unexpected data size.');
              assert(N == size(d.Y, 2), 'Unexpected output size.');
              assert(d.funId == func, 'Unexpected function id.');
              assert(inst == 0 || d.instId == inst, 'Unexpected instance id.');
              data{instance} = d;
            end

            % model types loop
            for modelType_cell = opts.modelTypes
              modelType = modelType_cell{:};

              % array of option designs to test for the current model type
              modelOpts = modelOptions.(modelType);
              modelOptsInd = modelOptionsInd.(modelType);

              % options loop
              for optInd = 1:length(modelOptsInd)
                modelOpt = modelOpts{optInd};
                modelOptInd = modelOptsInd(optInd);

                % prepare output files
                % hash = modelHash(modelOpts);
                res_fname = sprintf(opts.res_fname_template, dim, func, inst, N, ...
                  designType, modelType, modelOptInd);
                res_dir = fullfile(opts.exppath, opts.res_path, sprintf('%dD', dim), sprintf('f%d', func));

                if ~exist(res_dir, 'dir')
                  mkdir(res_dir);
                end

                res_fname = fullfile(res_dir, res_fname);

                if exist(res_fname, 'file')
                  if opts.rewrite_results
                    warning('Overwriting file ''%s''', res_fname);
                  else
                    warning('File ''%s'' exists, skipping the CV loop.', ...
                      res_fname);
                    continue;
                  end
                end

                % sampling random parameters
                if isfield(opts.rndModelOptions, modelType)
                  rng('shuffle');
                  paramDef = opts.rndModelOptions.(modelType);

                  assert(isfield(opts.rndFullFact, modelType));
                  fullFact = opts.rndFullFact.(modelType);

                  k = randi(size(fullFact, 1));
                  j = 1;
                  for param_cell = fields(paramDef)'
                    param = param_cell{:};
                    vals = paramDef.(param);
                    if iscell(vals)
                      modelOpt.(param) = vals{fullFact(k, j)};
                    else
                      modelOpt.(param) = vals(fullFact(k, j));
                    end
                    j = j + 1;
                  end

                  if opts.use_rng_seed
                    rng(opts.rng_seed);
                  end
                else
                  k = NaN;
                end

                try
                  results = testOneModel(data, dim, func, inst_samp, N, ...
                    modelType, modelOpt, ...
                    cv, opts.cv_ind);
                catch err
                  if ~strcmp(err.identifier, 'testOneModel:tr_failed')
                    rethrow(err);
                  else
                    continue;
                  end
                end

                if ~isempty(results) && (opts.rewrite_results || ~exist(res_fname, 'file'))
                  sleep_time = 1;
                  for t = 1:opts.save_retries
                    pause(sleep_time);
                    iserr = 0;
                    try
                      save(res_fname, 'dim', 'func', 'inst', 'inst_samp', 'N_cell', 'N', 'designType', ...
                        'cv', 'results', 'modelType', 'modelOptInd', 'modelOpt', 'k');
                    catch err
                      iserr = 1;
                      report = getReport(err);
                      warning('save failed with error:\n%s', report);
                    end

                    if ~iserr, break; end
                    sleep_time = 1.5 * sleep_time;
                  end % save retries

                  if iserr, error('Could not save file within %d retries.', opts.save_retries); end
                end
              end % options loop

            end % model types loop
          end % design types loop
        end % dataset sizes loop

      end % instance loop
    end % function loop
  end % dimension loop

  % clean up parpool
%  if opts.use_parpool
%    MetaParPool('close');
%  end

end % function


function results = testOneModel(data, dim, func, inst_samp, N, ...
  modelType, modelOptions, ...
  cv, cv_ind)

  fprintf('======= Cross-validation =======\n');
  fprintf('      # of folds:  %d\n', length(cv_ind));
  fprintf('      function:    %d\n', func);
  fprintf('      dimension:   %d\n', dim);
  fprintf('      instances:   %s\n', num2str(inst_samp));
  fprintf('      data size:   %d\n', N);
  fprintf('      model type:  %s\n', modelType);
  fprintf('      model opts:\n');
  disp(modelOptions);
  fprintf('================================\n');

  % the result structure, one row per each CV fold
  c = numel(cv_ind);
  res = struct('model', cell(c, 1), 'mse', cell(c, 1), 'mae', cell(c, 1), 'r2', cell(c, 1), 'n', cell(c, 1));

  % parallel loop over specified CV folds
  for i = cv_ind
    if numel(inst_samp) == 1
      d = data{inst_samp(1)};
    else
      assert(i <= numel(inst_samp));
      d = data{inst_samp(i)};
    end

    assert(~isempty(d));
    assert(isstruct(d) && isfield(d, 'X') && isfield(d, 'Y'));

    xmean = zeros(1, dim);
    generation = 0;

    Xtr = d.X(:, cv.training(i));
    Ytr = d.Y(:, cv.training(i));
    Xte = d.X(:, cv.test(i));
    Yte = d.Y(:, cv.test(i));

    try
      mdl = ModelFactory.createModel(modelType, modelOptions, xmean, generation);
      mdl = mdl.trainModel(Xtr', Ytr', xmean, generation);

      % the following might cause result files to be large
      % res(i).model = mdl;

      Ypred = mdl.modelPredict(Xte');
      n = numel(Ypred);

      se = norm(Yte - Ypred, 2).^2
      mse = se / n;
      mae = norm(Yte - Ypred, 1) / n;
      res(i).mse = mse;
      res(i).mae = mae;
      res(i).r2 = 1 - se / var(Yte);
      res(i).n = n;
      res(i).inst = inst_samp(i);
      res(i).Ypred = Ypred;
    catch err
      report = getReport(err);
      warning(['Training of model ''%s'' on %dD, func %d, inst %s, N %d' ...
               ' failed with error:\n%s'], ...
        modelType, dim, func, num2str(inst_samp), N, report);
      err = MException('testOneModel:tr_failed', ['Training of ' num2str(i) 'th fold failed']);
      throw(err);
    end
  end % CV loop

  results = res;

end % function


function cv = make_cvpartition(opts, N)
  % rng init
  if opts.use_rng_seed
    rng(opts.rng_seed);
  end

  cv = cvpartition(N, opts.cv_type, opts.cv_param);

  % check that cv indices specified in options are valid
  assert(all(1 <= opts.cv_ind & opts.cv_ind <= cv.NumTestSets), ...
    sprintf('Invalid cross-validation indices: %s.', num2str(opts.cv_ind)));
end


function res = myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
