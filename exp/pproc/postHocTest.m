function [pv, summary] = postHocTest(data, varargin)
%MULTCOMP Pairwise comparison of algorithms on a set of benchmarking
%   functions. Performs Friedman posthoc test and adjusts pairwise p-values
%   error by Bergmann-Hommel correction.
%   Currently implemented as a wrapper over an R script.
%
%   Input:
%     data - a table with benchmarking function as rows and algorithms as
%            variables containing best achieved fitness values
%     posthoc_test - 'Friedman' (default) | 'Quade' | 'FriedmanAlignedRanks'
%
%   Output:
%     pv      - a table of adjusted p-values with algorithms in columns and
%            rows
%     summary - test summary such as Friedman's mean ranks
%
%   TODO:
%     [ ] comparison of multiple algorithms with a control

  Rpath = fullfile('exp', 'pproc', 'R');
  Rscript = 'postHocTest.R';
  fin = [tempname, '.in.csv'];
  fout_pv = [tempname, '.pv.out.csv'];
  fout_sum = [tempname, '.sum.out.csv'];
  
  % invert the data so that higher is better
  u = max(max(data));
  l = min(min(data));
  data = (u - data) / (u - l);

  if nargin >= 2
    if ismember(varargin{1}, ...
      {'friedman', 'quade', 'aligned ranks'})
      test = varargin{1};
    else
      error('Unrecognized test ''%s''', varargin{1});
    end
  else
    test = 'friedman';
  end

  % write the input file with floats in scientific notation
  dlmwrite(fin, data, 'precision', '%e', 'delimiter', ',');

  args = sprintf(' -i "%s" -p "%s" -s "%s" --posthoc_test "%s"', ...
    fin, fout_pv, fout_sum, test);
  cmd = [fullfile(Rpath, Rscript), args];

  % execute the script and check exit code
  [status, cmdout] = system(cmd);

  if status ~= 0
    delete(fin);
    if exist(fout_pv, 'file'), delete(fout_pv), end
    if exist(fout_sum, 'file'), delete(fout_sum), end

    error('Script %s ended with error status %d and output:\n%s', ...
      Rscript, status, cmdout);
  end
  
  % read the script's output or outputs
  pv = dlmread(fout_pv, ',');
  summary = dlmread(fout_sum, ',');

  % clean up
  delete(fin);
  if exist(fout_pv, 'file'), delete(fout_pv), end
  if exist(fout_sum, 'file'), delete(fout_sum), end
end

