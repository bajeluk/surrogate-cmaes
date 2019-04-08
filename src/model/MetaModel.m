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
    metaTrainFlag         % flag denoting result of training selection
    archive
    population
    XTrain
    yTrain
    modelsCount           % number of models in the MetaModel
    models                % instances of models
    trainedModels         % indicatior of trained models
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

      obj.models = [];
      obj.isModelTrained = false;
      obj.dim       = size(xMean, 2);
      obj.shiftMean = zeros(1, obj.dim);
      obj.shiftY    = 0;
      obj.stdY      = 1;
      
      obj.dataset.X = [];
      obj.dataset.y = [];

      % general model prediction options
      obj.predictionType = defopts(modelOptions, 'predictionType', 'fValues');
      obj.transformCoordinates = defopts(modelOptions, 'transformCoordinates', true);
      obj.dimReduction = defopts(modelOptions, 'dimReduction', 1);
      obj.options.normalizeY = defopts(modelOptions, 'normalizeY', true);

      % meta model
      obj.models = cell(1, obj.modelsCount);
      obj.trainedModels = false(1, obj.modelsCount);
      
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
      if (~isempty(X) && ~isempty(y))
        obj.dataset.X = X;
        obj.dataset.y = y;
      end
      obj.population = population;
      
      obj.trainSigma = stateVariables.sigma;
      obj.trainBD = stateVariables.BD;
      
      % model training
      obj.models = cell(1, obj.modelsCount);
      obj.trainedModels = false(1, obj.modelsCount);
      
      % select best model id
      [bestModelId, obj.metaTrainFlag] = obj.chooseBestModelTrain();
      
      % train all selected models
      for m = 1:numel(bestModelId)
        % copy general and new values from the chosen settings
        metaModelOpts = obj.metaModelOptions;
        bestSetFields = fieldnames(obj.metaModelOptions.parameterSets(bestModelId(m)));
        for f = 1:numel(bestSetFields)
          metaModelOpts.(bestSetFields{f}) = ...
            obj.metaModelOptions.parameterSets(bestModelId(m)).(bestSetFields{f});
        end

        try
          % model constructor
          newModel = GpModel(metaModelOpts, obj.xMean);
          % train new model
          newModel = newModel.train(obj.getDataset_X(), obj.getDataset_y(), obj.stateVariables, ...
                                    obj.sampleOpts, obj.archive, obj.population);

          % save new model
          obj.models{bestModelId(m)} = newModel;
          obj.trainedModels(bestModelId(m)) = true;
        catch
          fprintf(2, 'MetaModel.train(): Selected model id %d training failed.\n', bestModelId(m));
        end
      end
            
      % if any model is trained, the metamodel can make predictions
      if any(obj.trainedModels)
        obj.isModelTrained = true;
      else
        obj.metaTrainFlag = -1;
      end
      % only one model trained
      if sum(obj.trainedModels) == 1
        obj.metaTrainFlag = 0;
      end
      
    end

    function [y, sd2] = modelPredict(obj, X)
    % MetaModel prediction method
      
      y = []; sd2 = [];
    
      % is another selection necessary?
      if sum(obj.trainedModels) > 1
        modelId = obj.chooseBestModelPredict(X);
      % only one model trained
      elseif sum(obj.trainedModels) == 1
        modelId = find(obj.trainedModels);
      % no model trained
      else
        fprintf(2, 'MetaModel.predict(): No model trained. MetaModel prediction failed.\n');
        return
      end
    
      % model prediction
      if obj.isTrained
        [y, sd2] = obj.models{modelId}.modelPredict(X);
      else
        fprintf(2, 'MetaModel.predict(): Model prediction failed.\n');
      end
      
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

    function [bestModelIndex, metaTrainFlag] = chooseBestModelTrain(obj)
    % Select the most convenient models for actual training dataset.
    
      switch lower(obj.bestModelSelection)
        case 'pitra2019landscape'
          [bestModelIndex, metaTrainFlag] = obj.getTree2019Train();
          % calculate metafeatures
          % find appropriate model
        otherwise
          error(['MetaModel.chooseBestModelTrain: ', ...
                 obj.MetaModelOptions.bestModelSelection, ...
                 ' -- no such option available']);
      end

    end
    
    function [bestModelIndex] = chooseBestModelPredict(obj, Xtest)
    % Select the most convenient model for actual prediction dataset.
    
      switch lower(obj.bestModelSelection)
        case 'pitra2019landscape'
          bestModelIndex = obj.getTree2019Predict(Xtest);
          % calculate metafeatures
          % find appropriate model
        otherwise
          error(['MetaModel.chooseBestModelPredict: ', ...
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

    function [modelId, metaTrainFlag] = getTree2019Train(obj)
    % Select the models using decision tree from article Pitra et al.
    % (2019): Landscape Analysis of Gaussian Process Surrogates for the
    % Covariance Matrix Adaptation Evolution Strategy
    
      % TODO: get proper bounds for decision splits
    
      % models used 1:LIN, 2:SE, 3:MATERN5, 4:RQ, 5:GIBBS
      modelId = []; 
      metaTrainFlag = 0; % only one model trained
      
      % data variables
      X_A = obj.archive.X;
      y_A = obj.archive.y;
      X_T = obj.getDataset_X();
      y_T = obj.getDataset_y();
      
      % feature settings
      settings.cma_cov = obj.trainBD*obj.trainBD';
      settings.cma_mean = obj.xMean;
      settings.cma_step_size = obj.trainSigma;
      
      % calculate CMA features on archive
      ft_cma_A = feature_cmaes(X_A, y_A, settings);
      % first node of the tree
      if ft_cma_A.cma_lik < -1.05034e8
        % necessary dispersion on traintest set will be calculated in 
        % getTree2019Predict
        modelId = [3, 4];
        metaTrainFlag = 1;
      else
        % calculate metamodel features on archive
        ft_mm_A = feature_ela_metamodel(X_A, y_A, settings);
        if ft_mm_A.lin_simple_adj_r2 < 0.965407
          if obj.dim < 15
            % calculate dispersion features on archive
            ft_dis_A = feature_dispersion(X_A, y_A, settings);
            if ft_dis_A.ratio_mean_25 < 0.245507
              % metamodel fts on archive
              if ft_mm_A.quad_simple_adj_r2 < 0.999398
                modelId = 3; % MATERN5
              else
                modelId = 4; % RQ
              end
            else
              % calculate metamodel fts on train set
              ft_mm_T = feature_ela_metamodel(X_T, y_T, settings);
              if ft_mm_T.quad_simple_adj_r2 < 0.998876
                % calculate levelset fts on archive
                ft_lvl_A = feature_ela_levelset(X_A, y_A, settings);
                if ft_lvl_A.mmce_qda_50 < 0.205074
                  modelId = 5; % GIBBS
                else
                  % metamodel fts on archive
                  if ft_mm_A.quad_simple_adj_r2 < 0.218322
                    % calculate dispersion fts on train set
                    ft_dis_T = feature_dispersion(X_T, y_T, settings);
                    if ft_dis_T.ratio_median_02 < 0.0313419
                      modelId = 2; % SE
                    else
                      modelId = 5; % GIBBS
                    end
                  else
                    % necessary information content fts on traintest set
                    % will be calculated in getTree2019Predict
                    modelId = [4, 5];
                    metaTrainFlag = 2;
                  end
                end
              else
                modelId = 4; % RQ
              end
            end
          else % dimension >= 15
            % calculate metamodel fts on train set
            ft_mm_T = feature_ela_metamodel(X_T, y_T, settings);
            if ft_mm_T.quad_simple_adj_r2 < 0.999377
              % metamodel fts on train set
              if ft_mm_T.quad_w_interact_adj_r2 < 1
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
      
    end % obj.getTree2019Train
    
    function modelId = getTree2019Predict(obj, Xtest)
    % Select the models using decision tree from article Pitra et al.
    % (2019): Landscape Analysis of Gaussian Process Surrogates for the
    % Covariance Matrix Adaptation Evolution Strategy
    
      X_Tp = [obj.getDataset_X(); Xtest];
      y_Tp = [obj.getDataset_y(); NaN(size(Xtest, 1), 1)];
    
      switch obj.metaTrainFlag
        case 1
          % calculate dispersion on traintest set
          ft_dis_Tp = feature_dispersion(X_Tp, y_Tp);
          if ft_dis_Tp.ratio_median_02 < 0.10666
            modelId = 3; % MATERN5
          else
            modelId = 4; % RQ
          end
        case 2
          % calculate information content fts on traintest set
          ft_inf_Tp = feature_infocontent(X_Tp, y_Tp);
          if ft_inf_Tp.eps_max < 0.0783446
            modelId = 5; % GIBBS
          else
            modelId = 4; % RQ
          end
        otherwise
          % return existing modelId's
          modelId = find(obj.trainedModels);
      end
    
    end
    
  end

  methods (Static)
    function result = calculateTrainRange(percentile, dimension)
      result = chi2inv(percentile,dimension);
    end
  end
end
