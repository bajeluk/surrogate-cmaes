classdef RegModel < Model
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
  end
    
    methods
        function obj = RegModel(modelOptions, xMean)
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
            mterms = defopts(modelOptions, 'modelterms', 1);
            if (numel(mterms) == 1)
                obj.modelterms = buildcompletemodel(mterms, numel(xMean));
            else
                obj.modelterms = mterms;
            end
            
        end
        
        function nData = getNTrainData(obj)
            % returns the required number of data for training the model
            nData = length(obj.modelterms) + 1;
        end
        
        function obj = trainModel(obj, X, y, xMean, generation, population, cmaesVariables)
            if (~isempty(X) && ~isempty(y))
                obj.dataset.X = X;
                obj.dataset.y = y;
            end
            obj.trainLikelihood = 0.1;
      
            obj.polyModel = polyfitn(X, y, obj.modelterms);
            obj.trainGeneration = generation;
        end

        function [y, sd2] = modelPredict(obj, X)
            y = polyvaln(obj.polyModel, X);
            sd2 = var(y);
        end
    end
end

