function [pv] = multComp(data)
%MULTCOMP Pairwise comparison of algorithms on a set of benchmarking
%   functions. Performs Friedman test and adjusts pairwise p-values
%   error by Bergmann-Hommel correction.
%   Currently implemented as a wrapper over an R script.
%
%   Input:
%     data - a table with benchmarking function as rows and algorithms as
%            variables containing best achieved fitness values
%
%   Output:
%     pv   - a table of corrected p-values with algorithms in columns and
%            rows
%
%   TODO:
%     [ ] comparison of multiple algorithms with a control

  Rpath = fullfile('exp', 'pproc', 'R');
  Rscript = 'multcomp.R';
  input = [tempname, '.csv'];
  output = [tempname, '.csv'];

  args = sprintf(' --input %s --output %s', input, output);
  cmd = [fullfile(Rpath, Rscript), args];

  % execute the script
  [status, cmdout] = system(cmd);

  if status ~= 0
    delete(input, cmdout);
    error('Script %s ended with error status %d and output:\n%s', ...
      Rscript, status, cmdout);
  end

  pv = readtable(output, ...
    'ReadVariableNames', true, ...
    'ReadRowNames', true, ...
    'Delimiter', ',' ...
  );

  delete(input, output);
end

