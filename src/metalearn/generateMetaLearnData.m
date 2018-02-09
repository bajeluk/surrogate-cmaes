function generateMetaLearnData(dims, funIds, instIds, Ns, designs, smooth, ...
    output_path, fname_template, rewrite_files)
%GENERATEMETALEARNDATA Generate data sets for meta learning.
%   A separate data set is created for each combination of the following
%   parameters:
%     dims           -- array of dimensionalities
%     funIds         -- array ids of BBOB functions
%     instIds        -- array of BBOB instance ids
%     Ns             -- cell array of expressions (strings) for data set sizes
%     designs        -- cell array of design types (strings)
%
%   Other inputs:
%     smooth         -- flag for smoothing designs (boolean)
%     output_path    -- the output folder
%     fname_template -- the template for the data set file name
%     rewrite_files  -- whether to rewrite existing files

  if ~exist(output_path, 'dir')
    mkdir(output_path)
  end

  % dimension loop
  for dim = dims
    dim_dir = fullfile(output_path, sprintf('%dD', dim));
    if ~exist(dim_dir, 'dir')
      mkdir(dim_dir);
    end

    % number of points loop
    for N_cell = Ns
      % design of points generating loop
      for design_cell = designs
        design = design_cell{:};
        N = myeval(N_cell{:});
        X = makeDesign(N, dim, smooth, design)';

        % function loop
        for funId = funIds
          % instances loop
          for instId = instIds
            Y = evalDesign(X, funId, instId);

            % save data
            fname = sprintf(fname_template, dim, funId, instId, N, design);
            fname = fullfile(dim_dir, fname);

            if exist(fname, 'file') && ~rewrite_files
              warning('Dataset already exists, skipping. File: %s', fname);
              continue;
            elseif exist(fname, 'file')
              warning('Rewriting file %s', fname);
            else
              fprintf('Saving file %s\n', fname);
            end

            save(fname, 'X', 'Y', 'funId', 'instId');
          end % instIds
        end % funIds
      end % designs
    end % Ns
  end % dims
end % function