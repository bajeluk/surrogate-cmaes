classdef LmmModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was built
    trainSigma            % sigma of the generation when the model was built
    trainBD               % BD of the generation when the model was built
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling
    options

    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space

    polyModel
    trainLikelihood
    modelterms

    kernel                % model kernel function
    knn                   % number of nearest neighbours
    modelTerms
    archive
    M
  end
    
  methods
    function obj = LmmModel(modelOptions, xMean)
      % constructor
      assert(size(xMean, 1) == 1, 'LmmModel (constructor): xMean is not a row-vector.');
      obj.options = modelOptions;

      % computed values
      obj.useShift  = defopts(obj.options, 'useShift', false);
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
      
      % local meta model settings
      obj.kernel = defopts(modelOptions, 'kernel', @(x) (1 - x.^2).^2);

       % Turn off warnings
%             warning('off','MATLAB:rankDeficientMatrix')
    end
        
    function nData = getNTrainData(obj)
      % returns the required number of data for training the model
      nData = obj.dim * (obj.dim + 3)/2 + 3;
    end
        
        function obj = trainModel(obj, X, y, xMean, generation, archive)
            if (~isempty(X) && ~isempty(y))
                obj.dataset.X = X;
                obj.dataset.y = y;
            end
            obj.trainLikelihood = 0.1;
            obj.trainGeneration = generation;
            
            
            % All we have to do is save archive, since lmm model cannot
            % really be trained
            [pointsInArchive, ~] = size(archive.X);
            obj.archive = archive;
            obj.knn = ceil(min(2*(getNTrainData(obj)-2), sqrt(pointsInArchive * (getNTrainData(obj) - 2))));
          
        end

    function [y, sd2] = modelPredict(obj, X)
      % predict function values for matrix of input points X using local
      % meta model (locally-weighted regression)

      % calculate M matrix from sigma and BD
      obj.M = calculateM(obj, obj.stateVariables);

      [nPoints, ~] = size(X);

      y = NaN(nPoints, 1);
      % query point loop
      for i = 1:nPoints
        queryPoint = X(i, :);
        % get points closest to the query point
        [~, closestY, closestZ, closestDistance] = getClosestPoints(obj, queryPoint);
        % predict function value for the query point
        y(i) = predictQueryPoint(obj, closestY, closestZ, closestDistance);
      end

      % set identical sd2 for all points (could be estimated)
      sd2 = var(y) * ones(nPoints, 1);
    end
        
        function [x] = minimumX(obj, archive)
            ub = obj.sampleOpts.ubounds;
            lb = obj.sampleOpts.lbounds;
        
            cmaesopt.LBounds = lb;
            cmaesopt.UBounds = ub;
            cmaesopt.SaveVariables = false;
            cmaesopt.LogModulo = 0;
            cmaesopt.DispModulo = 0;
            cmaesopt.DispFinal = 0;
            cmaesopt.Seed = 'inherit';
            sigma = [0.3*(ub - lb)];
            % sigma(end) = min(10*mean(sigma(1:end-1)), sigma(end));
            % there is ARD covariance
            % try run cmaes for 500 funevals to get bounds for covariances
            cmaesopt.MaxFunEvals = 500;
        
            eval_func = @(X) obj.predict(X');
            [opt, fval] = s_cmaes(eval_func, obj.trainMean, sigma, cmaesopt);        
            x = opt';
        end
  end

  methods (Access = protected)
    function [closestX, closestY, closestZ, closestDistance] = getClosestPoints(obj, queryPoint)
    % get points from the Archive closest to the query point

      % get Z values of X points (change of variables)
      zValues = (obj.M * (obj.archive.X - repmat(queryPoint, size(obj.archive.X, 1), 1))')';
      % Mahalanobis distance is Euclidean distance with the new
      % variable Z computing Euclidean distance
      distances = sqrt(sum(zValues.*zValues, 2));
      % find the closest points by sorting distance
      [~, I] = sort(distances);
      % select only k-nearest points from the archive, their Z values
      % and distances
      selectedIndexes = I(1:obj.knn);
      closestX = obj.archive.X(selectedIndexes, :);
      closestY = obj.archive.y(selectedIndexes, :);
      closestZ = zValues(selectedIndexes, :);
      closestDistance = distances(selectedIndexes);
    end
        
    function prediction = predictQueryPoint(obj, y, zValues, distances)
    % predict query point fitness value

      % feature length for full quadratic model
      lsqDim = obj.dim*(obj.dim + 3)/2 + 1;

      % compose Ztilda from z values in dimension D
      % Ztilda =
      %   (z_1^2, ..., z_D^2, z_1*z_2, ..., z_{D-1}*z_D, z_1, ..., z_D, 1)
      % using 'ones' function instantly determinates constant term to 1
      Ztilda = ones(obj.knn, lsqDim);
      jj=1;
      % quadratic terms
      Ztilda(:, jj:jj+obj.dim-1) = zValues.^2;
      jj = jj + obj.dim;
      % interaction terms
      for j = 2:obj.dim
        Ztilda(:, jj:jj+j-2) = repmat(2 * zValues(:, j), 1, j-1) .* zValues(:, 1:j-1);
        jj = jj + j - 1;
      end
      % linear terms
      Ztilda(:, jj:jj+obj.dim-1) = zValues;

      % calculate weights
      W = sqrt(obj.kernel(distances / distances(obj.knn)));

      WZ = repmat(W, 1, lsqDim) .* Ztilda;
      Wy = W .* y;

      % Beta = (W*Z)\(W*y);
      Beta = WZ \ Wy;
      % estimate of the function value at the query point corresponds
      % to the last coefficient of Beta, i.e., constant term of the
      % quadratic model
      prediction = Beta(end);
    end
        
        function M = calculateM(obj, cmaesVariables)
            D = diag(cmaesVariables.diagD);
            B = cmaesVariables.BD / D;
            
            
            
            M = diag(1./sqrt(diag(obj.stateVariables.sigma^2*D.^2)))*B'; %'
            
            %Dpow = D^(-1/2);
            %Btrans = transpose(B);
            %M = (1 / cmaesVariables.sigma) * Dpow * Btrans;
        end
  end
end

