classdef MetaModel < Model
  properties    % derived from abstract class "Model"
    dim                   % dimension of the input space X (determined from x_mean)
    trainGeneration = -1; % # of the generation when the model was built
    trainMean             % mean of the generation when the model was trained
    trainSigma            % sigma of the generation when the model was trained
    trainBD               % BD of the generation when the model was trained
    dataset               % .X and .y
    useShift = false;
    shiftMean             % vector of the shift in the X-space
    shiftY = 0;           % shift in the f-space
    predictionType        % type of prediction (f-values, PoI, EI)
    transformCoordinates  % transform X-space
    stateVariables        % variables needed for sampling new points as CMA-ES do
    sampleOpts            % options and settings for the CMA-ES sampling

    % GpModel specific fields
    stdY                  % standard deviation of Y in training set, for normalizing output
    options
    hyp
    meanFcn
    covFcn
    likFcn
    infFcn
    nErrors
    trainLikelihood

    % Dimensionality-reduction specific fields
    dimReduction          % Reduce dimensionality for model by eigenvectors
    % of covatiance matrix in percentage
    reductionMatrix       % Matrix used for dimensionality reduction

    % MetaModel specific properties
    metaModelOptions
    archive
    population
    XTrain
    yTrain
    modelsCount           % number of models in the MetaModel
    model                 % instances of models, 2D cell array of modelscount*historyLength+1
    bestModelIndex
    bestModelSelection    % method of best model selection
    choosingCriterium     % values of calculated criterium (mse/mae/...) for each model
    nTrainData            % min of getNTrainData of all created models
    xMean
    isModelTrained
  end

  methods (Access = public)
    function obj = MetaModel(modelOptions, xMean)
    % MetaModel constructor
    
      obj.metaModelOptions = modelOptions;
      obj.xMean = xMean;
      obj.modelsCount = length(modelOptions.parameterSets);
      assert(obj.modelsCount ~= 0, 'MetaModel(): No model provided!');

      obj.bestModelSelection = defopts(modelOptions, 'bestModelSelection', 'pitra2019landscape');

      obj.model = [];
      obj.isModelTrained = false;
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      obj.stdY      = 1;

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.dimReduction = defopts(modelOptions, 'dimReduction', 1);
      obj.options.normalizeY = defopts(modelOptions, 'normalizeY', true);

      obj.nTrainData = 3*obj.dim;
    end

    function nData = getNTrainData(obj)
      nData = obj.nTrainData;
    end

    function trained = isTrained(obj)
      % check whether the model chosen as the best in the newest generation is trained
      if (isempty(obj.isModelTrained))
        trained = false;
      else
        trained = obj.isModelTrained;
      end
    end

    function obj = trainModel(obj, ~, ~, ~, ~)
      % This function is empty because it is not needed, training is
      % done in train().
    end

    function obj = train(obj, X, y, stateVariables, sampleOpts, archive, population)
    % New train function for MetaModel
    
      obj.archive = archive;
      obj.stateVariables = stateVariables;
      obj.sampleOpts = sampleOpts;
      obj.xMean = stateVariables.xmean';
      obj.XTrain = X;
      obj.yTrain = y;
      obj.population = population;
      
      obj.trainSigma = stateVariables.sigma;
      obj.trainBD = stateVariables.BD;
      
      obj.isModelTrained = true;
      
      
      
    end

    function [y, sd2] = modelPredict(obj, X)
    % MetaModel prediction method
    
      % select best model id
      bestModelId = obj.chooseBestModel(X);
      % copy general and new values from the chosen settings
      metaModelOpts = obj.metaModelOptions;
      bestSetFields = fieldnames(obj.metaModelOptions.parameterSets(bestModelId));
      for f = 1:numel(bestSetFields)
        metaModelOpts.(bestSetFields{f}) = ...
          obj.metaModelOptions.parameterSets(bestModelId).(bestSetFields{f});
      end
      
      % model constructor
      newModel = GpModel(metaModelOpts, obj.xMean);
      % train new model
      newModel = newModel.train(obj.XTrain, obj.yTrain, obj.stateVariables, ...
                                obj.sampleOpts, obj.archive, obj.population);
      if ~newModel.isTrained && ~isempty(obj.model)
        newModel = obj.model.train(obj.XTrain, obj.yTrain);
        if ~newModel.isTrained
          newModel = obj.model;
        end
      end
      
      % model prediction
      if newModel.isTrained
        obj.model = newModel;
      
        [y, sd2] = obj.model.modelPredict(X);
      else
        y = []; sd2 = [];
        fprintf(2, 'MetaModel.predict(): Model training failed. No prediction available.');
      end
      
    end

    function X = getDataset_X(obj)
      X = obj.models{obj.bestModelIndex,1}.getDataset_X();
    end

    function y = getDataset_y(obj)
      y = obj.models{obj.bestModelIndex,1}.getDataset_y();
    end

  end

  methods (Access = protected)
    function gpModel = createGpModel(obj, modelIndex, xMean)
      if (isstruct(obj.MetaModelOptions.parameterSets))
        newModelOptions = obj.MetaModelOptions.parameterSets(modelIndex);
      elseif (iscell(obj.MetaModelOptions.parameterSets))
        newModelOptions = obj.MetaModelOptions.parameterSets{modelIndex};
      end
      newModelOptions.predictionType = obj.predictionType;
      newModelOptions.transformCoordinates = obj.transformCoordinates;
      newModelOptions.dimReduction = obj.dimReduction;
      newModelOptions.options.normalizeY = obj.options.normalizeY;

      gpModel = ModelFactory.createModel('gp', newModelOptions, xMean);

    end

    function [bestModelIndex] = chooseBestModel(obj, Xtest)
    % Select the most convenient model for actual training plus testing 
    % dataset.
    
      switch lower(obj.bestModelSelection)
        case 'pitra2019landscape'
          bestModelIndex = obj.getTree2019(Xtest);
          % calculate metafeatures
          % find appropriate model
        otherwise
          error(['MetaModel.chooseBestModel: ', ...
                 obj.MetaModelOptions.bestModelSelection, ...
                 ' -- no such option available']);
      end

    end

    function obj = copyPropertiesFromBestModel(obj)
      obj.stdY = obj.models{obj.bestModelIndex,1}.stdY;
      obj.options = obj.models{obj.bestModelIndex,1}.options;
      obj.hyp = obj.models{obj.bestModelIndex,1}.hyp;
      obj.meanFcn = obj.models{obj.bestModelIndex,1}.meanFcn;
      obj.covFcn = obj.models{obj.bestModelIndex,1}.covFcn;
      obj.likFcn = obj.models{obj.bestModelIndex,1}.likFcn;
      obj.infFcn = obj.models{obj.bestModelIndex,1}.infFcn;
      obj.nErrors = obj.models{obj.bestModelIndex,1}.nErrors;
      obj.trainLikelihood = obj.models{obj.bestModelIndex,1}.trainLikelihood;

      obj.shiftY = obj.models{obj.bestModelIndex,1}.shiftY;
      obj.trainSigma = obj.models{obj.bestModelIndex,1}.trainSigma;
      obj.trainBD = obj.models{obj.bestModelIndex,1}.trainBD;
      obj.trainMean = obj.models{obj.bestModelIndex,1}.trainMean;
      obj.shiftMean = obj.models{obj.bestModelIndex,1}.shiftMean;
      obj.reductionMatrix = obj.models{obj.bestModelIndex,1}.reductionMatrix;
      obj.stateVariables = obj.models{obj.bestModelIndex,1}.stateVariables;
      obj.sampleOpts = obj.models{obj.bestModelIndex,1}.sampleOpts;
    end

    function modelId = getTree2019(obj, Xtest)
    % Select the model using decision tree from article Pitra et al.
    % (2019): Landscape Analysis of Gaussian Process Surrogates for the
    % Covariance Matrix Adaptation Evolution Strategy
    
      % TODO: get proper bounds for decision splits
    
      % models used 1:LIN, 2:SE, 3:MATERN5, 4:RQ, 5:GIBBS
      modelId = []; 
      
      % data variables
      X_A = obj.archive.X;
      y_A = obj.archive.y;
      X_T = obj.XTrain;
      y_T = obj.yTrain;
      X_Tp = [X_T; Xtest];
      y_Tp = [y_T; NaN(size(Xtest, 1), 1)];
      
      % feature settings
      settings.cma_cov = obj.trainBD*obj.trainBD';
      settings.cma_mean = obj.xMean;
      settings.cma_step_size = obj.trainSigma;
      
      % calculate CMA features on archive
      ft_cma_A = feature_cmaes(X_A, y_A, settings);
      % first node of the tree
      if ft_cma_A.cma_lik < 0
        % calculate dispersion on traintest set
        ft_dis_Tp = feature_dispersion(X_Tp, y_Tp, settings);
        if ft_dis_Tp.ratio_median_02 < 0
          modelId = 3; % MATERN5
        else
          modelId = 4; % RQ
        end
      else
        % calculate metamodel features on archive
        ft_mm_A = feature_ela_metamodel(X_A, y_A, settings);
        if ft_mm_A.lin_simple_adj_r2 < 0
          if obj.dim < 15
            % calculate dispersion features on archive
            ft_dis_A = feature_dispersion(X_A, y_A, settings);
            if ft_dis_A.ratio_mean_25 < 0
              % metamodel fts on archive
              if ft_mm_A.quad_simple_adj_r2 < 0
                modelId = 3; % MATERN5
              else
                modelId = 4; % RQ
              end
            else
              % calculate metamodel fts on train set
              ft_mm_T = feature_ela_metamodel(X_T, y_T, settings);
              if ft_mm_T.quad_simple_adj_r2 < 0
                % calculate levelset fts on archive
                ft_lvl_A = feature_ela_levelset(X_A, y_A, settings);
                if ft_lvl_A.mmce_qda_50 < 0
                  modelId = 5; % GIBBS
                else
                  % metamodel fts on archive
                  if ft_mm_A.quad_simple_adj_r2 < 0
                    % calculate dispersion fts on train set
                    ft_dis_T = feature_dispersion(X_T, y_T, settings);
                    if ft_dis_T.ratio_median_02 < 0
                      modelId = 2; % SE
                    else
                      modelId = 5; % GIBBS
                    end
                  else
                    % calculate information content fts on traintest set
                    ft_inf_Tp = feature_infocontent(X_Tp, y_Tp, settings);
                    if ft_inf_Tp.eps_max < 0
                      modelId = 5; % GIBBS
                    else
                      modelId = 4; % RQ
                    end
                  end
                end
              else
                modelId = 4; % RQ
              end
            end
          else % dimension >= 15
            % calculate metamodel fts on train set
            ft_mm_T = feature_ela_metamodel(X_T, y_T, settings);
            if ft_mm_T.quad_simple_adj_r2 < 0
              % metamodel fts on train set
              if ft_mm_T.quad_w_interact_adj_r2 < 0
                modelId = 4; % RQ
              else
                modelId = 3; % MATERN5
              end
            else
              modelId = 4; % RQ
            end
          end
        else
          modelId = 1; % LIN
        end
      end
      
    end % obj.getTree2019
    
  end

  methods (Static)
    function result = calculateTrainRange(percentile, dimension)
      result = chi2inv(percentile,dimension);
    end
  end
end
