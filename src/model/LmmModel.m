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
            
            obj.modelTerms = buildcompletemodel(2, obj.dim);
            %Turn off warnings
            warning('off','MATLAB:rankDeficientMatrix')
        end
        
        function nData = getNTrainData(obj)
            % returns the required number of data for training the model)
            nData = obj.dim*(obj.dim+3)/2 + 2;
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
            obj.M = calculateM(obj, obj.stateVariables);
            obj.knn = ceil(min(2*getNTrainData(obj)-2, sqrt(pointsInArchive * (getNTrainData(obj) - 1 ))));
          
        end

        function [y, sd2] = modelPredict(obj, X)
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
                W = [W; sqrt(kernel(obj, (distances(i) / distances(obj.knn))))];
            end
            W = diag(W);
            
            [points, ~] = size(x);
            [terms, dim] = size(obj.modelTerms);
            Z = zeros(points, terms);
            for i = 1:terms
                for j = 1:dim
                    exponent = obj.modelTerms(i, j);
                    if exponent > 0
                        res = zValues(j, :).^exponent;
                        Z(:, i) = Z(:, i) + transpose(res);
                    end
                end
            end
            Z(:, terms) = 1;
            
            Beta = (W*Z)\(W*y);
            prediction = Beta(terms);
        end
        
        function M = calculateM(obj, cmaesVariables)
            D = diag(cmaesVariables.diagD);
            B = cmaesVariables.BD * inv(D);
            Dpow = D^(-1/2);
            Btrans = transpose(B);
            M = (1 / cmaesVariables.sigma) * Dpow * Btrans;
        end
    end
end

