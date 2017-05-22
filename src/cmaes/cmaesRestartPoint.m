function x0 = cmaesRestartPoint(archive, varargin)
% x0 = cmaesRestartPoint generates new restart point for the CMA-ES
% algorithm

  % initialize
  if nargin < 1
    x0 = [];
    return
  end
  restartSettings = settings2struct(varargin);
  
  dim = archive.dim;
  
  nSamples = defopts(restartSettings, 'NumSamples', 100);
  nSubsamplePoints = defopts(restartSettings, 'NumArchiveSamples', 5000);
  distRatio = defopts(restartSettings, 'DistanceRatio', 1/3);
  lb = defopts(restartSettings, 'LBound', zeros(1, dim));
  ub = defopts(restartSettings, 'UBound', ones(1, dim));
  
  % generate uniformly distributed points in sample space
  xRes = rand(nSamples, dim) .* repmat(ub - lb, nSamples, 1) + repmat(lb, nSamples, 1);
  % scatter(xRes(:,1), xRes(:,2), 'x', 'g')
  % hold on
  
  % load full archive data
  X = archive.X;
  nArchPoints = size(X, 1);
  % scatter(X(:,1), X(:,2), '+', 'b')
  
  % sample nSubsample archive points
  if nSubsamplePoints > nArchPoints
    nSubsamplePoints = nArchPoints;
  end
  archId = randsample(nArchPoints, nSubsamplePoints);
  XSub = X(archId, :);
  
  % calculate distances from newly sampled points to archive points and
  % return a fraction of points with the largest distance
  xResSubDist = sum(pdist2(xRes, XSub, 'euclidean'), 2);
  xNewId = find(xResSubDist >= quantile(xResSubDist, 1-distRatio));
  
  % return randomly chosen point with the large distance from the archive
  % points
  x0 = xRes(xNewId(randi(length(xNewId))), :);
  % scatter(x0(:,1), x0(:,2), 'o', 'r')
  % hold off

end