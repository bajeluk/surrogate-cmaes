function [pv, summary] = postHocTest(data, varargin)
%MULTCOMP Pairwise comparison of algorithms on a set of benchmarking
%   functions. Performs Friedman posthoc test and adjusts pairwise p-values
%   error by Bergmann-Hommel correction.
%   Currently implemented as a wrapper over an R script.
%
%   Input:
%     data - a table with benchmarking function as rows and algorithms as
%            variables containing best achieved fitness values
%     posthoc_test - 'friedman' (default) | 'quade' | 'aligned ranks'
%     (Friedman Aligned Ranks)
%     correction - 'bergmann' (default) | 'shaffer' | 'bonferroni'
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

  if nargin >= 2
    if ismember(varargin{1}, ...
      {'friedman', 'wilcoxon', 'quade', 'aligned ranks'})
      test = varargin{1};
    else
      error('Unrecognized test: ''%s''', varargin{1});
    end
  else
    test = 'friedman';
  end

  if nargin >= 3
    if ismember(varargin{2}, ...
      {'bergmann', 'shaffer', 'bonferroni'})
      corr = varargin{2};
    else
      error('Unrecognized correction type: ''%s''', varargin{2});
    end
  else
    corr = 'bergmann';
  end

  % write the input file with floats in scientific notation
  dlmwrite(fin, data, 'precision', '%e', 'delimiter', ',');

  args = sprintf(' -i "%s" -p "%s" -s "%s" -t "%s" -c "%s"', ...
    fin, fout_pv, fout_sum, test, corr);
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

