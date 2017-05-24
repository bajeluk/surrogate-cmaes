function x0 = cmaesRestartPoint(archive, varargin)
% x0 = cmaesRestartPoint(archive) generates new restart point for the 
% CMA-ES algorithm using S-CMA-ES archive of originally evaluated points.
%
% x0 = cmaesRestartPoint(dim) generates new restart point uniformly from
% [0, 1].
%
% Input:
%   archive - instance of class Archive
%   settings - pairs of property (string) and value or struct with 
%              properties as fields:
%     
%     'DistanceRatio'     - fraction of most distant points used to 
%                           randomly select the restart point | [0, 1]
%     'LBound'            - lower bounds of restart search space | double
%     'NumArchiveSamples' - number of sampled archive points (set Inf to
%                           use all the archive points) | integer
%     'NumSamples'        - number of sampled restart points | integer
%     'PlotResult'        - plot all sampled points | boolean
%     'UBound'            - upper bounds of restart search space | double
%
% Output:
%   x0 - CMA-ES restart point
%
% See Also:
%   s_cmaes

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
  nSubsamplePoints = defopts(restartSettings, 'NumArchiveSamples', Inf);
  distRatio = defopts(restartSettings, 'DistanceRatio', 1/3);
  lb = defopts(restartSettings, 'LBound', zeros(dim, 1));
  ub = defopts(restartSettings, 'UBound', ones(dim, 1));
  plotResult = defopts(restartSettings, 'plotResult', false);
  
  if isnumeric(archive)
    x0 = rand(dim, 1) .* (ub - lb) + lb;
    return
  end
  
  % generate uniformly distributed points in sample space
  xRes = rand(dim, nSamples) .* repmat(ub - lb, 1, nSamples) + repmat(lb, 1, nSamples);
  
  % load full archive data
  X = archive.X;
  nArchPoints = size(X, 1);
  
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
  xResSubDist = min(pdist2(xRes', XSub, 'euclidean'), [], 2);
  xNewId = find(xResSubDist >= quantile(xResSubDist, 1-distRatio));
  
  % return randomly chosen point with the large distance from the archive
  % points
  x0 = xRes(:, xNewId(randi(length(xNewId))));
  
  % plot results of this function (mainly for debugging purposes)
  if plotResult
    % combinations of dimensions
    dimComb = combnk(1:dim, 2);
    nComb = size(dimComb, 1);
    for dc = 1:nComb
      figure(dc)
      scatter(xRes(dimComb(dc, 1), :), xRes(dimComb(dc, 2), :), 'x', 'g')
      hold on
      scatter(X(:, dimComb(dc, 1)), X(:, dimComb(dc, 2)), '+', 'b')
      scatter(x0(dimComb(dc, 1), :), x0(dimComb(dc, 2), :), 'o', 'r')
      hold off
    end
  end

end