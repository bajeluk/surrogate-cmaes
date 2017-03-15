function [pv] = multComp(data, varargin)
%MULTCOMP Pairwise comparison of algorithms on a set of benchmarking
%   functions. Performs Friedman test and adjusts pairwise p-values
%   error by Bergmann-Hommel correction.
%   Currently implemented as a wrapper over an R script.
%
%   Input:
%     data - a table with benchmarking function as rows and algorithms as
%            variables containing best achieved fitness values
%     test - 'Friedman' (default) | 'Quade' | 'FriedmanAlignedRanks'
%
%   Output:
%     pv   - a table of corrected p-values with algorithms in columns and
%            rows
%
%   TODO:
%     [ ] comparison of multiple algorithms with a control

  Rpath = fullfile('exp', 'pproc', 'R');
  Rscript = 'multcomp.R';
  fin = [tempname, '.in.csv'];
  fout = [tempname, '.out.csv'];

  if nargin >= 2
    if ismember(varargin{1}, ...
      {'Friedman', 'Quade', 'FriedmanAlignedRanks'})
      test = varargin{1};
    else
      error('Unrecognized test ''%s''', varargin{1});
    end
  else
    test = 'Friedman';
  end

  args = sprintf(' -i %s -o %s --posthoc_test %s', fin, fout, test);
  cmd = [fullfile(Rpath, Rscript), args];

  % write the input file with floats in scientific notation
  dlmwrite(fin, data, 'precision', '%e', 'delimiter', ',');

  % execute the script and check exit code
  [status, cmdout] = system(cmd);

  if status ~= 0
    delete(fin);
    if exist(fout, 'file'), delete(fout), end

    error('Script %s ended with error status %d and output:\n%s', ...
      Rscript, status, cmdout);
  end

  % read the script's output
  pv = dlmread(fout, ',');

  % clean up
  delete(fin);
  if exist(fout, 'file'), delete(fout), end
end

