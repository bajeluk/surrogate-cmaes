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
        
        knn
        modelTerms
        archive
        M
    end
    
    methods
        function obj = LmmModel(modelOptions, xMean)
            % constructor
            assert(size(xMean,1) == 1, 'RegModel (constructor): xMean is not a row-vector.');
            obj.options = modelOptions;
      
            % computed values
            obj.useShift  = defopts(obj.options, 'useShift', false);
            obj.dim       = size(xMean, 2);
            obj.shiftMean = zeros(1, obj.dim);
            obj.shiftY    = 0;
      
            % general model prediction options
            obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
            obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
            
            %Turn off warnings
            warning('off','MATLAB:rankDeficientMatrix')
        end
        
        function nData = getNTrainData(obj)
            % returns the required number of data for training the model)
            nData = obj.dim*(obj.dim+3)/2 + 3;
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
            %sigma^2*D.^2,B,options.kernel,knn
            
            obj.M = calculateM(obj, obj.stateVariables);
            [points, dim] = size(X);
            y = [];
            for i = 1:points
                queryPoint = X(i, :);
                [closestX, closestY, closestZ, closestDistance] = getClosestPoints(obj, queryPoint);
                y = [y; predictQueryPoint(obj, closestX, closestY, closestZ, closestDistance)];
            end
            
            sd2 = var(y);
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
            zValues = calculateZValues(obj, obj.archive.X, queryPoint);
            distances = vecnorm(zValues);
            
            [sorted, I] = sort(distances);
            selectedIndexes = I(1:obj.knn);
            closestX = obj.archive.X(selectedIndexes, :);
            closestY = obj.archive.y(selectedIndexes, :);
            closestZ = zValues(:, selectedIndexes);
            closestDistance = distances(selectedIndexes);
        end
        
        function zValues = calculateZValues(obj, archiveX, queryPoint)
            [len, dim] = size(archiveX);
            zValues = [];
            for i = 1:len
                zValues = [zValues zValue(obj, archiveX(i,:), queryPoint)];
            end
        end
        
        function val = zValue(obj, p1, p2)
            val = obj.M * transpose(p1 - p2);
        end
       
        function funcValue = kernel(obj, x)
            funcValue = (1 - x^2)^2;
        end
        
        function prediction = predictQueryPoint(obj, x, y, zValues, distances)
            W = [];
            for i = 1:obj.knn
                W = [W; kernel(obj, (distances(i) / distances(obj.knn)))];
            end
            W = diag(W);
            
            [points, ~] = size(x);
            dim = obj.dim;
            
            zValues = zValues';
            Ztilda = ones(points, dim*(dim+3)/2+1);
            jj=1;
            Ztilda(:, jj:jj+dim-1) = zValues.^2;
            jj = jj + dim;

            for j=2:dim
                Ztilda(:, jj:jj+j-2) = repmat(2 * zValues(:, j), 1, j-1) .* zValues(:, 1:j-1);
                jj = jj + j - 1;
            end
    
            Ztilda(:, jj:jj+dim-1) = zValues;
            jj = jj + dim;

   

            
            
            
            
            %Z2 = ones(points, terms);
            %Z2(:, 1:dim) = zValues'.^2;
            %i = dim+1;
            %for j=2:dim
            %    Z2(:, i) = zValues(j, :) .* 2 .* zValues(j-1, :);
            %    i = i + 1;
            %end
            %Z2(:, i:i+dim-1) = zValues';

            wz = sqrt(diag(W)) .* Ztilda;
            wy = sqrt(diag(W)) .* y;
            
            Beta = wz \ wy;
            %Beta = (W*Z)\(W*y);
            prediction = Beta(dim*(dim+3)/2+1);
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

