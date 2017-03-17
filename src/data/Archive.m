classdef Archive < handle
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
      if (~isempty(obj.X))
        isNotYetSaved = ~(ismember(X, obj.X, 'rows'));
        % TODO: put here some TolX criterion, not this silly 'ismember()'
      else
        isNotYetSaved = true(size(X,1),1);
      end
      obj.X = [obj.X; X(isNotYetSaved,:)];
      obj.y = [obj.y; y(isNotYetSaved,:)];
      if (length(generation) == 1)
        obj.gens = [obj.gens; generation * ones(sum(isNotYetSaved),1)];
      else
        generation = generation(:);     % ensure column vector
        obj.gens = [obj.gens; generation(isNotYetSaved)];
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
          randp = randperm(length(closerThan2SigmaIdx));
          X = [X; obj.X(closerThan2SigmaIdx(randp(1:n)),:)];
          y = [y; closerDataY(randp(1:n))];
        end
      end
    end

    function [X, y] = getClosestDataFromPoints(obj, n, xInput, sigma, BD, trainRange)
      % returns union of 'n'-tuples of points which are closest to each
      % of data points from the points in 'xInput'
      % using (sigma*BD)-metric
      % if (n == 0), all the available data are returned
      nData = length(obj.y);
      X = []; y = [];

      if (nData == 0)
        return;
      end

      if (~exist('trainRange','var'))
        trainRange=Inf;
      end

      if (nData > n)
        % there are more data than 'n'
        indicesToReturn = false(size(obj.y,1),1);

        % compute coordinates in the (sigma*BD)-basis
        BDinv = inv(sigma*BD);
        % for each point from xInput:
        for i = 1:size(xInput,1)
          xTransf = ( BDinv * (obj.X - repmat(xInput(i,:),nData,1))' )';
          diff2 = sum(xTransf.^2, 2);
          isInRange = diff2 < (trainRange ^ 2);
          % take up to 'n' closest points from current point xInput(i,:)
          [~, closest] = sort(diff2);
          closest((n+1):end) = [];

          % union these points with the previous points, if they are in
          % trainRange
          indicesToReturn(closest) = true;
        end
      else
        indicesToReturn = true(size(obj.y,1),1);
      end

      % return the final points
      X = obj.X(indicesToReturn,:);
      y = obj.y(indicesToReturn);
    end

    function [X, y] = getClosestDataFromPopulation(obj, n, population, trainRange, sigma, BD)
      % returns points which are closest to each
      % of data points from the points in 'population'
      % using (sigma*BD)-metric
      % if (n == 0), all the available data are returned
      nData = length(obj.y);
      X = []; y = [];

      if (nData == 0)
        return;
      end

      if (nData > n)
        % there are more data than 'n'
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
        isInRange = diff2 < (trainRange ^ 2);
        [~, closest] = sort(diff2);
        closest((n+1):end) = [];
        % take only points that are in trainRange
        indicesToReturn(closest) = true;
      else
        indicesToReturn = true(size(obj.y,1),1);
      end

      % return the final points
      X = obj.X(indicesToReturn,:);
      y = obj.y(indicesToReturn);
    end

    function [X, y] = getTrainsetData(obj, trainsetType, trainsetSizeMax, xMean, trainRange, sigma, BD, population)

      % if the trainRange is specified as quantile from Chi2,
      % convert the number to the metric value
      if (isnumeric(trainRange) && trainRange >= 0.9 && trainRange <= 1.0)
        trainRange = sqrt(chi2inv(trainRange, obj.dim));
      end

      switch lower(trainsetType)
        case 'allpoints'
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
          [X, y] = obj.getClosestDataFromPoints(trainsetSizeMax, xMean, sigma, BD, trainRange);
        case 'nearesttopopulation'
          [X, y] = obj.getClosestDataFromPopulation(trainsetSizeMax, population, trainRange, sigma, BD);
      end
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
