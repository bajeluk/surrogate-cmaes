funIds = 1:2;
dims = [5];
instIds = [1:2];
Ns = {'50 * dim'};
design = {'lhs', 'lhsnorm'};

exppath = fullfile('exp', 'experiments');
input_path = fullfile(exppath, 'data_metalearn');
in_fname_template = strjoin({'data_metalearn_', ...
  '%dD_', ...
  'f%d_', ...
  'inst%d_', ...
  'N%d_', ...
  'design-%s.mat'}, '');

output_path = fullfile(exppath, 'data_metalearn_fts');
out_fname_template = 'metafeatures_N-%s_design-%s.mat';

t0 = tic;
for N_cell = Ns
  for design_cell = design
    des = design_cell{:};
    % 3d cell for results; N and design type will be distinguished by file name
    mfts = cell(max(dims), max(funIds), max(instIds));

    for dim = dims
      for funId = funIds
        for instId = instIds
          % debug
          fprintf('%dD, f%d, inst%d ...\n', dim, funId, instId);

          % load input data
          N = myeval(N_cell{:});
          in_fname = sprintf(in_fname_template, dim, funId, instId, N, des);
          in_fname = fullfile(input_path, sprintf('%dD', dim), in_fname);
          data = load(in_fname);

          % compute metafeatures
          opts.lb = -5 * ones(1, dim);
          opts.ub = 5 * ones(1, dim);
          opts.features = {'basic', 'cm_angle', 'cm_convexity', ...
                   'cm_gradhomo', 'dispersion', 'ela_distribution', ...
                   'ela_levelset', 'ela_metamodel', 'infocontent', ...
                   'nearest_better', 'pca'};
          [res.ft, res.values] = getMetaFeatures(data.X', data.Y', opts);
          mfts{dim, funId, instId} = res;

          % debug
          fprintf('Elapsed time: %.2f sec.\n', (tic - t0) / 1e6);
        end
      end
    end % dim loop

    % save results
    Nstr = strrep(N_cell{:}, ' * ', '');
    out_fname = sprintf(out_fname_template, Nstr, des);
    out_fname = fullfile(output_path, out_fname);
    save(out_fname, 'funIds', 'dims', 'instIds', 'Ns', 'design', 'mfts');
  end
end
