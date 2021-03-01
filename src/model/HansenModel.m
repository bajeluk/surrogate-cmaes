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

    modelTerms            % individual terms of model polynomial
    modelParams           % parameteres (values) of polynomial model terms
    type                  % type of model polynomial (linear, quadratic, full-quadratic)

    trainLikelihood       % TODO: delete this property
  end
    
  methods
        function obj = HansenModel(modelOptions, xMean)
            % constructor
            assert(size(xMean,1) == 1, 'HansenModel (constructor): xMean is not a row-vector.');
            obj.options = modelOptions;
            
            obj.options.maxScalingFactor = defopts(obj.options, 'maxScalingFactor', 1);
      
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
            nData = 3;
        end
        
    function obj = trainModel(obj, X, y, ~, generation)
    % train the linear-quadratic model based on the data (X,y)

      % save training data
      if (~isempty(X) && ~isempty(y))
        obj.dataset.X = X;
        obj.dataset.y = y;
      end
      % TODO: delete trainLikelihood property (now neccessary due to stats)
      obj.trainLikelihood = NaN;
      obj.trainGeneration = generation;
      % create terms of polynomial model according to the number of points
      nPoints = size(X, 1);
      obj = createModelTerms(obj, nPoints);
      % get polynomial model values for recent X
      Z = modelValues(obj, X);

      weights = linspace(obj.options.maxScalingFactor, 1, nPoints);

      scaledY = weights' .* y;

      scaledZ = Z .* repmat(weights', 1, size(Z, 2));

      obj.modelParams = pinv(scaledZ) * scaledY;
    end
        
    function [obj] = createModelTerms(obj, nPoints)
    % create terms of polynomial model according to the number of points

      obj.modelTerms = [];
      % linear terms
      obj.type = 'linear';
      for i = 1 : obj.dim + 1
        vals = zeros(obj.dim + 1, 1);
        vals(i)= 1;
        obj.modelTerms = [obj.modelTerms, vals];
      end
      % (pure) quadratic terms
      if (nPoints > (obj.dim * 2 + 1) * 1.1)
        obj.type = 'quad';
        for i = 2 : obj.dim + 1
          vals = zeros(obj.dim + 1, 1);
          vals(i) = 2;
          obj.modelTerms = [obj.modelTerms, vals];
        end
      end
      % interaction terms (full-quadratic)
      if (nPoints > ((obj.dim * 2 + 1) + ((obj.dim / 2) * (obj.dim - 1))) * 1.1)
        obj.type = 'full';
        for i = 2 : obj.dim + 1
          for j = i + 1 : obj.dim + 1
            vals = zeros(obj.dim + 1, 1);
            vals(i) = 1;
            vals(j) = 1;
            obj.modelTerms = [obj.modelTerms, vals];
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
                            if (x(i, j)) == 0
                                x(i, j) = 1;
                            end
                            x(i, j) = x(i, j) * X(i, k).^obj.modelTerms(k, j);
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
        
    function x = minimumX(obj)
    % get input values of model minimum

      % linear polynom
      if strcmp(obj.type, 'linear')
        x = obj.dataset.X(end, :) - 2 * transpose(obj.modelParams(2:end));
      % quadratic polynom
      else
        k = 2 * obj.dim + 2;
        hessian = zeros(obj.dim, obj.dim);
        for i = 1 : obj.dim
          hessian(i, i) = obj.modelParams(i + obj.dim + 1);
          % full-quadratic polynom
          if strcmp(obj.type, 'full')
            for j = i + 1 : obj.dim
              hessian(i, j) = obj.modelParams(k) / 2;
              hessian(j, i) = obj.modelParams(k) / 2;
              k = k + 1;
            end
          end
        end
        baseCoeff = obj.modelParams(2:obj.dim + 1);
        x = (pinv(hessian) * (baseCoeff / -2))';
      end
    end
  end
end
