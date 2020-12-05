classdef HansenModel < Model
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
    
    modelTerms
    modelParams
    powers
    
    trainLikelihood
  end
    
    methods
        function obj = HansenModel(modelOptions, xMean)
            % constructor
            assert(size(xMean,1) == 1, 'HansenModel (constructor): xMean is not a row-vector.');
            obj.options = modelOptions;
      
            % computed values
            obj.useShift  = defopts(obj.options, 'useShift', false);
            obj.dim       = size(xMean, 2);
            obj.shiftMean = zeros(1, obj.dim);
            obj.shiftY    = 0;
      
            % general model prediction options
            obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
            obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', false);
            
            %mterms = defopts(modelOptions, 'modelterms', 1);
            %if (numel(mterms) == 1)
            %    obj.modelterms = buildcompletemodel(mterms, numel(xMean));
            %else
            %    obj.modelterms = mterms;
            %end
            
        end
        
        function nData = getNTrainData(obj)
            % returns the required number of data for training the model
            nData = 1;
        end
        
        function obj = trainModel(obj, X, y, xMean, generation, population, cmaesVariables)
            if (~isempty(X) && ~isempty(y))
                obj.dataset.X = X;
                obj.dataset.y = y;
            end
            obj.trainLikelihood = 0.1;
            obj.trainGeneration = generation;
            obj.modelTerms = createModelTerms(obj, X);
            
            Z = modelValues(obj, X);
            
            obj.modelParams = pinv(Z) * y;
            
        end
        
        function [modelTerms] = createModelTerms(obj, X)
            [points, ~] = size(X);
            modelTerms = [];
            for i=1:obj.dim+1
                vals = zeros(obj.dim + 1, 1);
                vals(i)= 1;
                modelTerms = [modelTerms vals];
            end
           
            if (points > (obj.dim * 2 + 1) * 1.1)
                for i=2:obj.dim+1
                    vals = zeros(obj.dim + 1, 1);
                    vals(i) = 2;
                    modelTerms = [modelTerms vals];
                end
            end
            if (points > ((obj.dim * 2 + 1) + ((obj.dim / 2) * (obj.dim - 1))) * 1.1)
                for i=2:obj.dim+1
                    for j=i+1:obj.dim+1
                        vals = zeros(obj.dim + 1, 1);
                        vals(i) = 1;
                        vals(j) = 1;
                        modelTerms = [modelTerms vals];
                    end
                end
            end
        end
        
        function x = modelValues(obj, X)
            [points, dim] = size(X);
            [rows, cols]= size(obj.modelTerms);
            X = [ones(points, 1) X];
            x = zeros(points, cols);
            for i=1:points
                for j=1:cols
                    for k=1:rows
                        if obj.modelTerms(k, j) ~= 0
                            x(i, j) = x(i, j) + X(i, k).^obj.modelTerms(k, j);
                        end
                    end
                end
            end
        end

        function [y, sd2] = modelPredict(obj, X)
            values = modelValues(obj, X);
            y = transpose(transpose(obj.modelParams) * transpose(values));
            sd2 = var(y);
        end
    end
end

