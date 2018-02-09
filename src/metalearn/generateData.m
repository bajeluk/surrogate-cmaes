rng(42); % for reproducibility

funIds = 1:24;
dims = [2, 5, 10, 20];
instIds = [1:5 40:50];
Ns = {'50 * dim'};
designs = {'lhs', 'lhsnorm'};
smooth = false;

exppath = fullfile('exp', 'experiments');
output_path = fullfile(exppath, 'data_metalearn');

fname_template = strjoin({'data_metalearn_', ...
  '%dD_', ...
  'f%d_', ...
  'inst%d_', ...
  'N%d_', ...
  'design-%s.mat'}, '');
rewrite_files = true;

generateMetaLearnData(dims, funIds, instIds, Ns, designs, ...
  smooth, output_path, fname_template, rewrite_files);