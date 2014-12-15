classdef SurrogateSelector
  methods (Static)
    function [xPreSample, zPreSample] = chooseDistantPoints(n, xSample, zSample, xTrain, xmean, sigma, BD)
      % choose 'n' of the xSample data, preferably not near to the xTrain data
      xPreSample = []; zPreSample = [];
      
      % compute coordinates in the (sigma*BD)-basis
      BDinv = inv(sigma*BD);
      xTrainTrans = ( BDinv * (xTrain - repmat(xmean',size(xTrain,1),1))' )';
      
      % cluster the xSample into 'n' clusters
      [idx, C] = kmeans(zSample, n);
      D = pdist2(C, xTrainTrans, 'euclidean');
      % find the nearest clusters from the points in 'xTrain'
      [~, nearestClusters] = min(D, [], 1);
      isXDistant = true(1,size(zSample,1));
      % omit the data of these clusters from xSample
      for cl = nearestClusters
        isXDistant(idx == cl) = 0;
      end
      % cluster the remaining points into 'n' clusters
      [idx2, C2, ~, D] = kmeans(zSample(isXDistant,:), n);
      % return the points nerest to the last 'n' centroids
      [~, closestToCentroid] = min(D, [], 1);
      xDistantIdx = find(isXDistant);
      xPreSample = xSample(xDistantIdx(closestToCentroid),:);
      zPreSample = zSample(xDistantIdx(closestToCentroid),:);
    end
  end
end
