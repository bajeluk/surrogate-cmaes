function algEvals = evaluateYThresholds(y_evals, thresholds)
% evaluateYThresholds -- calculate # function evaluations to get to y-values
%
% load the specified file, goes through the y_evals
% and calculate the # evaluations needed to reach f-values
% given by   10 .^ thresholds
% if threshold is not reached, NaN is returned in these entries

  % if (nargin >= 3)
  %   path = varargin{1};
  %   filename = fullfile(path, filename);
  % end
  %
  % load(filename, 'y_evals');

  % thresholds = [0 -1 -2 -3 -4 -5 -6 -7 -8];

  algEvals = nan(length(thresholds), length(y_evals));

  % for each algorithm runs (y_evals' fields)
  for i = 1:length(y_evals)
    % index into thresholds vector
    th = 1;
    % go through the y_evals matrix, where 
    %   1st column: best-f, 
    %   2nd column: # evaluations
    for ye = 1:length(y_evals{i})
      % go through the thresholds while the actual
      % so-far best-f is lower then the threshold
      while (th <= length(thresholds) ...
             && y_evals{i}(ye, 1) < 10^thresholds(th) ...
             && isnan(algEvals(th,i)))
        % save the # evaluations to resulting algEvals
        algEvals(th, i) = y_evals{i}(ye, 2);
        th = th + 1;
      end
    end
  end
end
