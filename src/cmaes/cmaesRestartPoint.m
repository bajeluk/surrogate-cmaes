function x0 = cmaesRestartPoint(archive, varargin)
% x0 = cmaesRestartPoint(archive) generates new restart point for the 
% CMA-ES algorithm using S-CMA-ES archive of originally evaluated points.
%
% x0 = cmaesRestartPoint(dim) generates new restart point uniformly from
% [0, 1].

  % initialize
  if nargin < 1
    x0 = [];
    return
  end
  restartSettings = settings2struct(varargin);
  
  if isnumeric(archive)
    dim = archive;
  else
    dim = archive.dim;
  end
  
  % parse settings
  nSamples = defopts(restartSettings, 'NumSamples', 100);
  nSubsamplePoints = defopts(restartSettings, 'NumArchiveSamples', 5000);
  distRatio = defopts(restartSettings, 'DistanceRatio', 1/3);
  lb = defopts(restartSettings, 'LBound', zeros(dim, 1));
  ub = defopts(restartSettings, 'UBound', ones(dim, 1));
  
  if isnumeric(archive)
    x0 = rand(dim, 1) .* (ub - lb) + lb;
    return
  end
  
  % generate uniformly distributed points in sample space
  xRes = rand(dim, nSamples) .* repmat(ub - lb, 1, nSamples) + repmat(lb, 1, nSamples);
  % scatter(xRes(:,1), xRes(:,2), 'x', 'g')
  % hold on
  
  % load full archive data
  X = archive.X;
  nArchPoints = size(X, 1);
  % scatter(X(:,1), X(:,2), '+', 'b')
  
  % empty archive
  if nArchPoints == 0
    x0 = xRes(:, 1);
    return
  end
  
  % sample nSubsample archive points
  if nSubsamplePoints > nArchPoints
    nSubsamplePoints = nArchPoints;
  end
  archId = randsample(nArchPoints, nSubsamplePoints);
  XSub = X(archId, :);
  
  % calculate distances from newly sampled points to archive points and
  % return a fraction of points with the largest distance
  xResSubDist = sum(pdist2(xRes', XSub, 'euclidean'), 2);
  xNewId = find(xResSubDist >= quantile(xResSubDist, 1-distRatio));
  
  % return randomly chosen point with the large distance from the archive
  % points
  x0 = xRes(:, xNewId(randi(length(xNewId))));
  % scatter(x0(:,1), x0(:,2), 'o', 'r')
  % hold off

end