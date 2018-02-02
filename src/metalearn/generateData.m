rng(42); % for reproducibility

funIds = 1:24;
dims = [2, 5, 10, 20];
instIds = [1:5 40:50];
Ns = {'50 * dim'};
exppath = fullfile('exp', 'experiments');
output_path = fullfile(exppath, 'data_metalearn');
rewrite_files = false;
designs = {'lhs'}; % {'lhs', 'lhsnorm'}

if ~exist(output_path, 'dir')
  mkdir(output_path)
end

fname_template = strjoin({'data_metalearn_', ...
                          '%dD_', ...
                          'instId%d_', ...
                          'N%d_', ...
                          'funId%d_', ...
                          'design-%s'}, '');

for dim = dims
  for Nstr = Ns
    for design = designs
      N = myeval(Nstr{:});
      X = makeDesign(N, dim, design{:})';

      for funId = funIds
        for instId = instIds
          Y = evalDesign(X, funId, instId);
          
          % save data
          fname = sprintf(fname_template, dim, instId, N, funId, design{:});
          fname = fullfile(output_path, fname);
          
          if exist(fname, 'file') && ~rewrite_files
            warning('Dataset %s already exists, skipping.');
            continue;
          elseif exist(fname, 'file')
            warning('Rewriting dataset %s');
          end
          
          save(fname, 'X', 'Y', 'funId', 'instId');
        end % instIds
      end % funIds
    end % designs
  end % Ns
end % dims
      