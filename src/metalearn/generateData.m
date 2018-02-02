rng(42); % for reproducibility

funIds = 1:24;
dims = [2, 5, 10, 20];
instIds = [1:5 40:50];
Ns = {'50 * dim'};
exppath = fullfile('exp', 'experiments');
output_path = fullfile(exppath, 'data_metalearn');
rewrite_files = true;
designs = {'lhs', 'lhsnorm'};
smooth = false;

if ~exist(output_path, 'dir')
  mkdir(output_path)
end

fname_template = strjoin({'data_metalearn_', ...
                          '%dD_', ...
                          'funId%d_', ...
                          'instId%d_', ...
                          'N%d_', ...
                          'design-%s.mat'}, '');

for dim = dims
  for N_cell = Ns
    for design_cell = designs
      design = design_cell{:};
      N = myeval(N_cell{:});
      X = makeDesign(N, dim, smooth, design)';

      for funId = funIds
        for instId = instIds
          Y = evalDesign(X, funId, instId);
          
          % save data
          fname = sprintf(fname_template, dim, instId, N, funId, design);
          fname = fullfile(output_path, fname);
          
          if exist(fname, 'file') && ~rewrite_files
            warning('Dataset already exists, skipping. File: %s', fname);
            continue;
          elseif exist(fname, 'file')
            warning('Rewriting file %s', fname);
          end
          
          save(fname, 'X', 'Y', 'funId', 'instId');
        end % instIds
      end % funIds
    end % designs
  end % Ns
end % dims
      