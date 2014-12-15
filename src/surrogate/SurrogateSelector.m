classdef SurrogateSelector
  properties

  end

  methods
    function [xPreSample, zPreSample] = chooseDistantPoints(n, xSample, zSample, xTrain, xmean, sigma, BD)
      % choose 'n' of the xSample data, preferably not near to the xTrain data
      nData = size(xSample, 1);
      xPreSample = []; zPreSample = [];
      
      if (nData == 0)
        return;
      end

      % compute coordinates in the (sigma*BD)-basis
      BDinv = inv(sigma*BD);
      xTrainTrans = ( BDinv * (xTrain - repmat(xmean',nData,1))' )';
      
      % cluster the xSample into 'n' clusters
      [idx, C] = kmeans(zSample, n);
      D = pdist2(C, xTrainTrans, 'euclidean');
      % find the nearest clusters from the points in 'xTrain'
      [~, nearestClusters] = min(D, 1);
      xDistantIdx = ones(1,size(zSample,1));
      % omit the data of these clusters from xSample
      for cl = nearestClusters
        xDistantIdx(idx == cl) = 0;
      end
      % cluster the remaining points into 'n' clusters
      [idx2, C2, ~, D] = kmeans(zSample(xDistantIdx,:), n);
      % return the points nerest to the last 'n' centroids
      
    end
  end
end
