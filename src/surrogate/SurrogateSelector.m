classdef SurrogateSelector
  methods (Static)
    function [xPreSample, zPreSample, xDistant] = chooseDistantPoints(n, xSample, zSample, xTrain, xmean, sigma, BD)
      % choose 'n' of the xSample data, preferably not near to the xTrain data
      %
      % n       -- the number of points to return
      % xSample, zSample (n x dim) -- random sample according to CMA-ES, row-vectors
      % xTrain  -- chosen archive evaluated points
      % xmean   -- current CMA-ES xmean
      % sigma, BD -- from CMA-ES
      %
      % returns column-wise vectors

      if (size(zSample,1) <= n)
        warning('SurrogateSelector.chooseDistantPoints(): not enough data to choose from!');
        xPreSample = xSample;
        zPreSample = zSample;
        return;
      end
      % xPreSample = []; zPreSample = [];

      if (size(xTrain,1) > 0)
        % compute coordinates in the (sigma*BD)-basis
        xTrainTrans = ( (sigma*BD) \ (xTrain - repmat(xmean',size(xTrain,1),1))' )';

        % cluster the xSample into 'n' clusters
        [idx, C] = kmeans(zSample, n);
        D = pdist2(C, xTrainTrans, 'euclidean');
        % find the nearest clusters from the points in 'xTrain'
        [~, nearestClusters] = min(D, [], 1);
        nearestClusters = unique(nearestClusters);
        isXDistant = true(1,size(zSample,1));
        % omit the data of these clusters from xSample
        % but leave at least 'n' points in isXDistant
        i = 1; cl = nearestClusters(i);
        while (sum(isXDistant) > n  &&  i <= length(nearestClusters))
          cl = nearestClusters(i);
          isXDistant(idx == cl) = 0;
          i = i + 1;
        end
        if (sum(isXDistant) < n)
          % if we do not have enough points, return back the last cluster
          isXDistant(idx == cl) = 1;
        end
        xDistant = xSample(isXDistant, :);
        % cluster the remaining points into 'n' clusters      
      else
        % cluster all the data in zSample
        isXDistant = true(size(zSample,1),1);
      end 
      [~, ~, ~, D] = kmeans(zSample(isXDistant,:), n);
      % return the points nerest to the last 'n' centroids
      [~, closestToCentroid] = min(D, [], 1);
      xDistantIdx = find(isXDistant);
      xPreSample = xSample(xDistantIdx(closestToCentroid),:)';
      zPreSample = zSample(xDistantIdx(closestToCentroid),:)';
    end

    function [xToReeval, xToReevalValid, zToReeval] = choosePointsToReevaluate(...
        xExtend, xExtendValid, zExtend, yExtend, nBest, nCluster)
      % choose points from extended population for reevaluation
      % with the original fitness

      % choose the best points
      [~, yExtendSortIdx] = sort(yExtend);
      bestIdx = yExtendSortIdx(1:nBest);
      xToReeval = xExtend(:,bestIdx);
      xToReevalValid = xExtendValid(:,bestIdx);
      zToReeval = zExtend(:,bestIdx);
      xExtend(:,bestIdx) = [];
      xExtendValid(:,bestIdx) = [];
      zExtend(:,bestIdx) = [];
      yExtend(bestIdx) = [];
      % cluster the rest of the points
      if (nCluster > 0)
        idx = kmeans(zExtend', nCluster);
        for cl = 1:nCluster
          % put the best from each cluster into 'xToReeval', 'zToReeval'
          thisClusterIdx = find(idx == cl);
          [~, clusterBestIdx] = min(yExtend(thisClusterIdx));
          clusterBestIdx = thisClusterIdx(clusterBestIdx);
          xToReeval = [xToReeval, xExtend(:,clusterBestIdx)];
          xToReevalValid = [xToReevalValid, xExtendValid(:,clusterBestIdx)];
          zToReeval = [zToReeval, zExtend(:,clusterBestIdx)];
        end
      end
    end
  end
end
