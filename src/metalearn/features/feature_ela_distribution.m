function ft = feature_ela_distribution(~, y)
% ft = feature_ela_distribution(X, y) returns ELA y-distribution features
% for dataset [X, y].
%
% Features:
%   skewness        - skewness of y values (R-package e1071 v1.6-8 type 3)
%   kurtosis        - kurtosis of y values (R-package e1071 v1.6-8 type 3)
%   number_of_peaks - number of y-distribution peaks

  if nargin < 2
    help feature_ela_distribution
    if nargout > 0
      ft = struct();
    end
    return
  end
  
  n = numel(y);

  % skewness according to R-package e1071 v1.6-8 skewness type 3
  ft.skewness = ((n-1)/n)^(3/2) * skewness(y, 1);
  % kurtosis according to R-package e1071 v1.6-8 kurtosis type 3
  ft.kurtosis = (1-1/n)^2 * kurtosis(y, 1) - 3;
  % number of peaks
  ft.number_of_peaks = numberOfPeaks(y);
end

function numOfPeaks = numberOfPeaks(x)
% Estimate the number of peaks in x distribution

  modemass_treshold = 0.01;
  % estimate density
  % TODO: proper settings of ksdensity function
  [y_dens, xi_dens] = ksdensity(x);
  n = length(y_dens);
  id = 2 : (n - 1);
  % find the position of valleys within the estimated y-distribution
  valleyId = [1, find(y_dens(id) < y_dens(id - 1) & y_dens(id) < y_dens(id + 1)) + 1, n + 1];
  % calculate the mass represented by a potential peak
  modemass = zeros(1, numel(valleyId) - 1);
  for i = 1:numel(valleyId) - 1
    modemass(i) = mean(y_dens(valleyId(i):valleyId(i+1) - 1)) * ...
                 abs(xi_dens(valleyId(i)) - xi_dens(valleyId(i+1) - 1));
  end
  numOfPeaks = sum(modemass > modemass_treshold);

end