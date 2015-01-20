classdef Archive
  properties
    dim  = 1;           % dimension of the input space X (determined from x_mean)
    X    = [];          % archive - input-space data
    y    = [];          % archive - dependent-space data
    gens = [];          % archive - generations of the data where they come from
  end

  methods
    function obj = Archive(dimension)
      % constructor
      obj.dim = dimension;
    end

    function obj = save(obj, X, y, generation)
      % save data (@X,@y) from @generation to the archive
      assert(size(X,1) == size(y,1), 'Archive.save: dimensions X and y mismatch.');
      assert(size(y,2) == 1, 'Archive.save: y is not a column vector.');
      assert(size(X,2) == obj.dim, 'Archive.save: dimension of X''s and Archive.dim mismatch.');
      obj.X = [obj.X; X];
      obj.y = [obj.y; y];
      obj.gens = [obj.gens; generation * ones(length(y),1)];
    end

    function [X, y] = getDataFromGenerations(obj, generations)
      % return data from generation(s) defined in scalar/vector @generations
      dataIdxs = ismember(obj.gens, generations);
      X = obj.X(dataIdxs, :);
      y = obj.y(dataIdxs);
    end

    function [X, y] = getDataNearPoint(obj, n, x, rangeSigma, sigma, BD)
      % returns up to 'n' data within distance of 'rangeSigma' along the point 'x'
      % using (sigma*BD)-metric
      % if more than 'n' data are closer than 'rangeSigma', k-means clustering is
      % performed
      % if (n == 0), all the available data are returned
      nData = length(obj.y);
      X = []; y = [];
      
      if (nData == 0)
        return;
      end

      % compute coordinates in the (sigma*BD)-basis
      BDinv = inv(sigma*BD);
      xTransf = ( BDinv * (obj.X - repmat(x,nData,1))' )';
      
      % take the points closer than *rangeSigma*
      diff = sum(xTransf.^2, 2);
      isInRange = diff < (rangeSigma ^ 2);

      if (sum(isInRange) <= n  ||  n <= 0)
        X = obj.X(isInRange,:);
        y = obj.y(isInRange);
      else
        % cluster the transformed data into n clusters
        closerDataX = xTransf(isInRange,:);
        closerDataY = obj.y(isInRange);
        closerThan2SigmaIdx = find(isInRange);
        [~, ~, ~, D] = kmeans(closerDataX, n);
        % D = ('n' x 'k') distances to the clusters' centroids
        % find the points nearest to the clusters' centers
        [~, closestToCentroid] = min(D, [], 1);
        for closestIdx = closestToCentroid
          % return the original coordinates, not the transformed
          X = [X; obj.X(closerThan2SigmaIdx(closestIdx),:)];
          y = [y; closerDataY(closestIdx)];
        end
      end
    end
  end
end
