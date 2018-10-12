function ft = feature_ela_distribution(~, y, settings)
% ft = FEATURE_ELA_DISTRIBUTION(X, y, settings) returns ELA y-distribution 
% features for dataset [X, y].
%
% The y-distribution features compute the kurtosis, skewness and number of 
% peaks of the kernel-based estimation of the density of the initial 
% design’s objective values. (Mersmann et al., 2011)
%
% settings:
%   modemass_treshold - treshold for the mass represented by a potential 
%                       peak | default: 0.01
%   skewness_type - 1: MATLAB skewness (flag = 1)
%                   2: MATLAB skewness (flag = 0)
%                   3: skewness according to R-package e1071 v1.6-8 
%                      skewness type 3 (default)
%   kurtosis_type - 1: MATLAB kurtosis (flag = 1)
%                   2: MATLAB kurtosis (flag = 0)
%                   3: skewness according to R-package e1071 v1.6-8 
%                      skewness type 3 (default)
%
% Features:
%   skewness        - skewness of y values
%   kurtosis        - kurtosis of y values
%   number_of_peaks - number of y-distribution peaks

  if nargin < 3
    if nargin < 2
      help feature_ela_distribution
      if nargout > 0
        ft = struct();
      end
      return
    end
    settings = struct();
  end
  
  % parse settings
  modemass_treshold = defopts(settings, 'modemass_treshold', 0.01);
  skewness_type = defopts(settings, 'skewness_type', 3);
  kurtosis_type = defopts(settings, 'kurtosis_type', 3);
  
  n = numel(y);

  % calculate skewness
  switch skewness_type
    case 1
      ft.skewness = skewness(y, 1);
    case 2
      ft.skewness = skewness(y, 0);
    case 3
      % skewness according to R-package e1071 v1.6-8 skewness type 3
      ft.skewness = ((n-1)/n)^(3/2) * skewness(y, 1);
    otherwise
      error('Wrong skewness_type')
  end
  
  % calculate kurtosis
  switch kurtosis_type
    case 1
      ft.kurtosis = kurtosis(y, 1);
    case 2
      ft.kurtosis = kurtosis(y, 0);
    case 3
      % kurtosis according to R-package e1071 v1.6-8 kurtosis type 3
      ft.kurtosis = (1-1/n)^2 * kurtosis(y, 1) - 3;
    otherwise
      error('Wrong kurtosis_type')
  end

  % number of peaks
  if all(isnan(y))
    ft.number_of_peaks = NaN;
  else
    ft.number_of_peaks = numberOfPeaks(y, modemass_treshold);
  end
end

function numOfPeaks = numberOfPeaks(x, modemass_treshold)
% Estimate the number of peaks in x distribution

  % estimate density, where the number of points is set to 512 according to 
  % the default value of density function in R and the method for bandwith
  % selection is probably different because there is "SJ" method used in
  % flacco (methods of Sheather & Jones (1991) to select the bandwidth 
  % using pilot estimation of derivatives)
  nPoints = 512;
  if verLessThan('matlab', 'R2017b')
    [y_dens, xi_dens] = ksdensity(x, 'npoints', nPoints);
  else
    [y_dens, xi_dens] = ksdensity(x, 'NumPoints', nPoints);
  end
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