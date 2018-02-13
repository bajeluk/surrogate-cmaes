function testMetaLearn(modelOptions, opts, funcs, dims, ...
  insts, Ns, models, designs)
%TESTMETALEARN Test models on metalearning data sets.

  if ~iscell(models), models = {models}; end
  if ~iscell(designs), designs = {designs}; end

  % Default options
  opts.cv_type = defopts(opts, 'cv_type', 'KFold');
  opts.cv_param = defopts(opts, 'cv_param', 10);
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

  % rng init
  if opts.use_rng_seed
    rng(opts.rng_seed);
  end

  % parpool init
  if opts.use_parpool
    poolobj = parpool(opts.parpool_size);
  else
    poolobj = [];
  end

  % CV partitioning
  cv = cvpartition(N, opts.cv_type, opts.cv_param);
  nTestSets = cv.NumTestSets;

  % dimension loop
  for dim = dims

    % function loop
    for func = funcs

      % instance loop
      for inst = insts

        % dataset sizes loop
        for N_cell = Ns

          % design types loop
          for design_cell = designs
            design = design_cell{:};

            % load the dataset
            fname = sprintf(opts.fname_template, dim, fun, inst, N, design);
            fname = fullfile(opts.dataset_path, sprintf('%dD', dim), fname);
            data = load(fname);

            % data sanity checks
            assert(dim == size(data.X, 1), 'Unexpected dimensionality.');
            assert(N == size(data.X, 2), 'Unexpected data size.');
            assert(N == size(data.Y, 2), 'Unexpected output size.');
            assert(data.funId == func, 'Unexpected function id.');
            assert(data.instId == inst, 'Unexpected instance id.');

            % models loop
            for model_cell = models
              model = model_cell{:};

              % options loop
              for opt_cell = modelOptions
                opt = opt_cell{:};
                testOneModel(dim, func, inst, N_cell, opt, 
              end % options loop

            end % models loop
          end % design loop
        end % dataset sizes loop
      end % instance loop
    end % function loop
  end % dimension loop

  % clean up parpool
  if ~isempty(poolobj)
    delete(poolobj)
  end

end % function


function testOneModel(dim, func, inst, N_cell, model, ...
)
  N = myeval(N_cell{:});

  % the result structure, one cell per CV fold
  cv_trained = cell(1, nTestSets);

  % prepare output files
  hash = modelHash(opt);
  res_fname = sprintf(opts.res_fname_template, dim, func, inst, N, ...
    design, model, hash);
  res_dir = fullfile(opts.exppath_short, opts.exp_id);

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
    end

    continue;
  end

  % parallel loop over CV folds
  parfor i = 1:nTestSets
    xmean = zeros(1, dim);
    generation = 0;

    Xtr = data.X(:, cv.train(i));
    Ytr = data.Y(:, cv.train(i));

    try
      mdl = ModelFactory.createModel(model, opt, xmean, generation);
      mdl = mdl.trainModel(Xtr', Ytr', xmean, generation);
      cv_trained{i} = mdl;
    catch err
      report = getReport(err);
      warning(['Training of model ''%s'' on %dD, func %d, inst %d,' ...
               ' failed with error:\n%s'], ...
        model, dim, func, inst, report);
    end
  end % CV loop

  if opts.rewrite_results || ~exist(res_fname, 'file')
    save(res_fname, 'dim', 'func', 'inst', 'N_cell', 'N', 'design', ...
      'X', 'Y', 'cv', 'cv_trained', 'model', 'opt');
  end
end % function


function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
