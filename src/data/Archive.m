classdef Archive < handle
  properties
    dim  = 1;           % dimension of the input space X (determined from x_mean)
    X    = [];          % archive - input-space data
    y    = [];          % archive - dependent-space data
    gens = [];          % archive - generations of the data where they come from
    tolX;               % minimal tolerance in X-space to consider points equal
                        % in maximal metric
  end

  methods
    function obj = Archive(dimension, varargin)
      % constructor
      obj.dim = dimension;
      if (nargin >= 2 && isnumeric(varargin{1}))
        obj.tolX = varargin{1};
      else
       obj.tolX = eps;
      end
    end

    function obj = save(obj, X, y, generation)
      % save data (@X,@y) from @generation to the archive
      assert(size(X,1) == size(y,1), 'Archive.save: dimensions X and y mismatch.');
      assert(size(y,2) == 1, 'Archive.save: y is not a column vector.');
      assert(size(X,2) == obj.dim, 'Archive.save: dimension of X''s and Archive.dim mismatch.');
      isNotYetSaved = ~obj.isInArchive(X);
      obj.X = [obj.X; X(isNotYetSaved,:)];
      obj.y = [obj.y; y(isNotYetSaved,:)];
      if (length(generation) == 1)
        obj.gens = [obj.gens; generation * ones(sum(isNotYetSaved),1)];
      else
        generation = generation(:);     % ensure column vector
        obj.gens = [obj.gens; generation(isNotYetSaved)];
      end
    end
    
    function obj = sortLast(obj, count)
        [points, ~] = size(obj.y);
        breakPoint = points - count;
        if breakPoint < 0
            breakPoint = 0;
        end
        toSort = obj.y(breakPoint + 1:points);
        [sortedVales, sortedIdxs] = sort(toSort, 'descend');
        
        newX = obj.X(1:breakPoint, :);
        newY = obj.y(1:breakPoint);
        
        for idx = sortedIdxs
            newX = [newX; obj.X(breakPoint + idx, :)];
            newY = [newY; obj.y(breakPoint + idx)];
        end 
        
        obj.X = newX;
        obj.y = newY;
    end

    function obj = delete(obj, indices)
      % remove data at given indices
      obj.X(indices,:) = [];
      obj.y(indices)   = [];
      obj.gens(indices) = [];
    end

    function [inArx, indices] = isInArchive(obj, X)
      inArx = false(size(X,1), 1);
      indices = zeros(size(X,1), 1);
      if (~isempty(obj.X))
        for i = 1:size(X, 1)
          [minDist, idx] = min(max(abs(repmat(X(i,:),size(obj.X,1),1) - obj.X), [], 2));
          inArx(i) = minDist <= obj.tolX;
          if (inArx(i))
            indices(i) = idx;
          end
        end
      end
    end

    function [X, y] = getDataFromGenerations(obj, generations)
      % return data from generation(s) defined in scalar/vector @generations
      dataIdxs = ismember(obj.gens, generations);
      X = obj.X(dataIdxs, :);
      y = obj.y(dataIdxs);
    end

    function [X, y, nData] = getDataNearPoint(obj, n, x, rangeSigma, sigma, BD)
      % returns up to 'n' data within distance of 'rangeSigma' along the point 'x'
      % using (sigma*BD)-metric
      % if more than 'n' data are closer than 'rangeSigma', k-means clustering is
      % performed
      % if (n == 0), all the available data are returned
      % returns:
      %   nData -- the number of all available data in the specified range

      MAX_POINTS_FOR_KMEANS = 1000;

      nData = length(obj.y);
      X = []; y = [];

      if (nData == 0)
        return;
      end

      % compute coordinates in the (sigma*BD)-basis
      xTransf = ( (sigma * BD) \ (obj.X - repmat(x,nData,1))' )';

      % take the points closer than *rangeSigma*
      diff2 = sum(xTransf.^2, 2);
      isInRange = diff2 < (rangeSigma ^ 2);
      nData = sum(isInRange);

      if (nData <= n  ||  n <= 0)
        X = obj.X(isInRange,:);
        y = obj.y(isInRange);
      else
        % cluster the transformed data into n clusters
        closerDataX = xTransf(isInRange,:);
        closerDataY = obj.y(isInRange);
        closerThan2SigmaIdx = find(isInRange);
        if (length(closerThan2SigmaIdx) <= MAX_POINTS_FOR_KMEANS)
          % Use k-means clustering if not too much points to choose from
          useClustering = true;
          try
            [~, ~, ~, D] = kmeans(closerDataX, n);
            % D = ('n' x 'k') distances to the clusters' centroids
            % find the points nearest to the clusters' centers
            [~, closestToCentroid] = min(D, [], 1);
            for closestIdx = closestToCentroid
              % return the original coordinates, not the transformed
              X = [X; obj.X(closerThan2SigmaIdx(closestIdx),:)];
              y = [y; closerDataY(closestIdx)];
            end
          catch err
            warning('Archive.getDataNearPoint(): %s\n', err.message);
            useClustering = false;
          end
        else
          useClustering = false;
        end
        if (~useClustering)
          randp = randperm(length(closerThan2SigmaIdx));
          X = [X; obj.X(closerThan2SigmaIdx(randp(1:n)),:)];
          y = [y; closerDataY(randp(1:n))];
        end
      end
    end

    function [X, y] = getClosestDataFromPoints(obj, n, xInput, sigma, BD, trainRange)
      % returns union of points which are closest to each
      % of data points from the points in 'xInput'
      % using (sigma*BD)-metric
      % at least  floor(n/size(xInput,1))  closest points are returned
      % from each point in xInput, maybe more
      %
      % if (n == 0), all the available data are returned
      nData = length(obj.y);
      X = []; y = [];

      if (nData == 0)
        return;
      end

      if (~exist('trainRange','var'))
        trainRange=Inf;
      end

      lambda = size(xInput,1);

      indicesToReturn = false(nData, 1);

      distanceOrderMatrix = zeros(lambda, nData);
      % compute coordinates in the (sigma*BD)-basis
      BDinv = inv(sigma*BD);

      % for each point from xInput:
      for i = 1:lambda
        xTransf = ( BDinv * (obj.X - repmat(xInput(i,:),nData,1))' )';
        diff2 = sum(xTransf.^2, 2);
        diff2Rank = ranking(diff2);
        isInRange = diff2 < (trainRange ^ 2);
        % fill the ranking of distances to the archive points from this xInput
        distanceOrderMatrix(i,isInRange) = diff2Rank(isInRange)';
      end

      nInRange = sum(any(distanceOrderMatrix > 0));
      if (nInRange > n)
        % there are more data than 'n' in the specified trainRange
        
        % initialize 
        nNearestPerXInput = floor(n / lambda);
        indicesToReturn(any( ...
            distanceOrderMatrix > 0 & distanceOrderMatrix <= nNearestPerXInput, 1)) ...
            = true;
        % increase the number of nearest points per xInput points
        % until we get enough points (ie. 'n' points) chosen
        while (sum(indicesToReturn) < n)
          nNearestPerXInput = nNearestPerXInput + 1;
          newIndicesToReturn = indicesToReturn;
          newIndicesToReturn(any( ...
              distanceOrderMatrix > 0 & distanceOrderMatrix <= nNearestPerXInput, 1)) ...
              = true;
          if (sum(newIndicesToReturn) > n)
            % we have got more points than we need, so choose
            % the ones with the best fitness
            nRest = n - sum(indicesToReturn);
            % choose from the set-diff of the last and current choice
            chooseFrom = find(newIndicesToReturn & ~indicesToReturn);
            yDiff = obj.y(chooseFrom);
            [~, syDiff] = sort(yDiff);
            indicesToReturn(chooseFrom(syDiff(1:nRest))) = true;
            break;
          end
          indicesToReturn = newIndicesToReturn;
        end
      else
        indicesToReturn(any(distanceOrderMatrix > 0)) = true;
      end

      % return the final points
      X = obj.X(indicesToReturn,:);
      y = obj.y(indicesToReturn);
    end

    function [X, y] = getClosestDataFromPopulation(obj, n, population, trainRange, sigma, BD)
      % returns points form Archive which are closest to any of
      % the points from 'population' using (sigma*BD)-metric
      % if (n == 0), all the available data are returned
      nData = length(obj.y);
      X = []; y = [];

      if (nData == 0)
        return;
      end

      indicesToReturn = false(size(obj.y,1),1);

      % compute coordinates in the (sigma*BD)-basis
      BDinv = inv(sigma*BD);
      % for each point from population:
      for i = 1:size(population.x,2)
        xTransf = ( BDinv * (obj.X - repmat(population.x(:,i)',nData,1))' )';
        if i==1
          diff2 = sum(xTransf.^2, 2);
        else
          diff2 = min(diff2, sum(xTransf.^2, 2));
        end
      end
      % take up to 'n' closest points
      isInRange = find(diff2 < (trainRange ^ 2));

      if (length(isInRange) > n)
        % there are more data than 'n'
        [~, closest] = sort(diff2(isInRange));
        closest((n+1):end) = [];
        % take only points that are in trainRange
        indicesToReturn(isInRange(closest)) = true;
      else
        indicesToReturn(isInRange) = true;
      end

      % return the final points
      X = obj.X(indicesToReturn,:);
      y = obj.y(indicesToReturn);
    end

    function [X, y] = getTrainsetData(obj, trainsetType, trainsetSizeMax, xMean, trainRange, sigma, BD, population, modelOptions)

      % if the trainRange is specified as quantile from Chi2,
      % convert the number to the metric value
      if (isnumeric(trainRange) && trainRange >= 0.9 && trainRange <= 1.0)
        trainRange = sqrt(chi2inv(trainRange, obj.dim));
      % if the trainRange is specified as small number
      % convert the number to multiple of 0.99 quantile of sqrt(Chi2 distribution)
      elseif (isnumeric(trainRange) && trainRange > 1.0 && trainRange <= 4.0)
        trainRange = trainRange * sqrt(chi2inv(0.99, obj.dim));
      end

      switch lower(trainsetType)
        case { 'allpoints', 'recent' }
          %get all data in range (0 as first parameter)
          [X, y] = obj.getDataNearPoint(0, xMean, trainRange, sigma, BD);
          if size(X,1)> trainsetSizeMax
            X = X(end-trainsetSizeMax+1:end,:);
            y = y(end-trainsetSizeMax+1:end,:);
            %remove elements from the beginning
          end
        case 'clustering'
          [X, y] = obj.getDataNearPoint(trainsetSizeMax, xMean, trainRange, sigma, BD);
        case 'nearest'
          [X, y] = obj.getClosestDataFromPoints(trainsetSizeMax, population.x', sigma, BD, trainRange);
        case 'nearesttopopulation'
          [X, y] = obj.getClosestDataFromPopulation(trainsetSizeMax, population, trainRange, sigma, BD);
        case 'latest'
          [X, y] = obj.getLatestData(trainsetSizeMax);
          if modelOptions.latestTruncateRatio > 0 && modelOptions.latestTruncateRatio < 1
            [points, ~] = size(X);
            points = ceil(points * modelOptions.latestTruncateRatio);
            
            [sortedValues, sortedIndexes] = sort(y);
            X = X(sortedIndexes(1:points), :);
            y = y(sortedIndexes(1:points));
          end
         
      end
    end
    
    function [X, y] = getLatestData(obj, trainsetSizeMax)
        [archiveSize, ~] = size(obj.y);
        dataCnt = min(trainsetSizeMax, archiveSize);
        X = obj.X(end-dataCnt+1:end, :);
        y = obj.y(end-dataCnt+1:end);
    end

    function obj2 = duplicate(obj)
      obj2 = Archive(obj.dim);
      obj2.X = obj.X;
      obj2.y = obj.y;
      obj2.gens = obj.gens;
    end

    function obj = restrictToGenerations(obj, toGens)
      id_gens = ismember(obj.gens, toGens);
      obj.X = obj.X(id_gens, :);
      obj.y = obj.y(id_gens, :);
      obj.gens = obj.gens(id_gens);
    end
  end
end
