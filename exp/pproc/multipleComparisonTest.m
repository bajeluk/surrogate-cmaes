function [pv, stat] = multipleComparisonTest(data, varargin)
%MULTIPLECOMPARISONTEST Multiple comparison of algorithms with either
%   Friedman or Iman-Davenport (test) test.
%   Calls an R script and reports back the test's statistic and p-value.

  Rpath = fullfile('exp', 'pproc', 'R');
  Rscript = 'multipleComparisonTest.R';
  
  fin = [tempname, '.in.csv'];
  fout_pv = [tempname, '.pv.out.csv'];
  fout_stat = [tempname, '.stat.out.csv'];

  if nargin >= 2
    if ismember(varargin{1}, {'friedman', 'iman'})
      test = varargin{1};
    else
      error('Unrecognized test ''%s''', varargin{1});
    end
  else
    test = 'iman';
  end

  % write the input file with floats in scientific notation
  dlmwrite(fin, data, 'precision', '%e', 'delimiter', ',');

  args = sprintf(' -i "%s" -p "%s" -s "%s" --test "%s"', ...
    fin, fout_pv, fout_stat, test);
  cmd = [fullfile(Rpath, Rscript), args];

  % execute the script and check exit code
  [status, cmdout] = system(cmd);

  if status ~= 0
    delete(fin);
    if exist(fout_pv, 'file'), delete(fout_pv), end
    if exist(fout_stat, 'file'), delete(fout_stat), end

    error('Script %s ended with error status %d and output:\n%s', ...
      Rscript, status, cmdout);
  end

  % read the script's output or outputs
  pv = dlmread(fout_pv, ',');
  stat = dlmread(fout_stat, ',');

  % clean up
  delete(fin);
  if exist(fout_pv, 'file'), delete(fout_pv), end
  if exist(fout_stat, 'file'), delete(fout_stat), end
end

