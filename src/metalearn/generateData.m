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
                          'f%d_', ...
                          'inst%d_', ...
                          'N%d_', ...
                          'design-%s.mat'}, '');
% dimension loop
for dim = dims
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
          fname = fullfile(output_path, fname);
          
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
      