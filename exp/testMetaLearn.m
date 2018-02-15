function testMetaLearn(modelOptions, modelOptionsInd, opts, funcs, dims, ...
  insts, Ns, modelTypes, designTypes)
%TESTMETALEARN Test models on metalearning data sets.

  if ~iscell(Ns), Ns = {Ns}; end
  if ~iscell(modelTypes), modelTypes = {modelTypes}; end
  if ~iscell(designTypes), designTypes = {designTypes}; end

  assert(isnumeric(funcs), '''funcs'' has to be numeric');
  assert(isnumeric(dims), '''dims'' has to be numeric');
  assert(isnumeric(insts), '''insts'' has to be numeric');

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

  % parpool init
  if opts.use_parpool
    poolobj = parpool(opts.parpool_size);
  else
    poolobj = [];
  end

  % dimension loop
  for dim = dims
    % function loop
    for func = funcs
      % instance loop
      for inst = insts

        % dataset sizes loop
        for N_cell = Ns
          N = myeval(N_cell{:});

          % CV partitioning
          cv = make_cvpartition(opts, N);

          % design types loop
          for design_cell = designTypes
            designType = design_cell{:};

            % load the dataset
            fname = sprintf(opts.fname_template, dim, fun, inst, N, designType);
            fname = fullfile(opts.dataset_path, sprintf('%dD', dim), fname);
            data = load(fname);

            % data sanity checks
            assert(dim == size(data.X, 1), 'Unexpected dimensionality.');
            assert(N == size(data.X, 2), 'Unexpected data size.');
            assert(N == size(data.Y, 2), 'Unexpected output size.');
            assert(data.funId == func, 'Unexpected function id.');
            assert(data.instId == inst, 'Unexpected instance id.');

            % model types loop
            for modelType_cell = modelTypes
              modelType = modelType_cell{:};

              % array of option designs to test for the current model type
              modelOpts = modelOptions.(modelType);
              modelOptsInd = modelOptionsInd.(modelType);

              % options loop
              for optInd = modelOptsInd
                modelOpt = modelOpts{optInd};

                % prepare output files
                % hash = modelHash(modelOpts);
                res_fname = sprintf(opts.res_fname_template, dim, func, inst, N, ...
                  designType, modelType, optInd);
                res_dir = fullfile(opts.exppath, 'results');

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

                results = testOneModel(dim, func, inst, N, modelType, modelOpt, cv);

                if ~isempty(results) && opts.rewrite_results || ~exist(res_fname, 'file')
                  save(res_fname, 'dim', 'func', 'inst', 'N_cell', 'N', 'designType', ...
                    'Y', 'cv', 'results', 'modelType', 'optInd', 'modelOpt');
                end
              end % options loop

            end % model types loop
          end % design types loop
        end % dataset sizes loop

      end % instance loop
    end % function loop
  end % dimension loop

  % clean up parpool
  if ~isempty(poolobj)
    delete(poolobj)
  end

end % function


function results = testOneModel(dim, func, inst, N, modelType, opt, cv)
  % the result structure, one row per each CV fold
  nTestSets = cv.NumTestSets;

  results = struct('model', [], 'Y', [], 'Ypred', []);

  % parallel loop over CV folds
  parfor i = 1:nTestSets
    xmean = zeros(1, dim);
    generation = 0;

    Xtr = data.X(:, cv.train(i));
    Ytr = data.Y(:, cv.train(i));
    Xte = data.X(:, cv.test(i));
    Yte = data.Y(:, cv.test(i));

    try
      mdl = ModelFactory.createModel(modelType, opt, xmean, generation);
      mdl = mdl.trainModel(Xtr', Ytr', xmean, generation);

      results(i).model = mdl;
      results(i).Y = Yte';
      results(i).Ypred = mdl.predict(Xte');
    catch err
      report = getReport(err);
      warning(['Training of model ''%s'' on %dD, func %d, inst %d, N %d' ...
               ' failed with error:\n%s'], ...
        model, dim, func, inst, N, report);
    end
  end % CV loop

end % function


function cv = make_cvpartition(opts, N)
  % rng init
  if opts.use_rng_seed
    rng(opts.rng_seed);
  end

  % CV partitioning
  cv = cvpartition(N, opts.cv_type, opts.cv_param);
end


function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
end
