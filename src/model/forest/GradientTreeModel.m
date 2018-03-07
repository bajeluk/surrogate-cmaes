classdef GradientTreeModel < TreeModel  
  properties %(Access = protected)
    tree_regularization % Ridge regularization parameter
  end
  
  methods
    function obj = GradientTreeModel(modelOptions)
      % constructor
      modelOptions.splitGain_weightedGain = false;
      modelOptions.tree_regularization = defopts(modelOptions, 'tree_regularization', 0);
      % splitGain_regularization and tree_regularization should have the
      % same values
      modelOptions.splitGain_regularization = defopts(modelOptions, ...
        'splitGain_regularization', modelOptions.tree_regularization);
      modelOptions.tree_splitGainFunc = @GradientSplitGain;
      modelOptions.tree_predictorFunc = @ConstantModel;
      modelOptions.tree_lossFunc = @mseLossFunc;
      obj = obj@TreeModel(modelOptions);
      
      % specific model options
      obj.tree_regularization = modelOptions.tree_regularization;
    end
  end
  
  methods (Access = protected)
    function predictor = trainPredictor(obj, X, y, modelID)
      G = sum(y(:, 1));
      H = sum(y(:, 2));
      w = repmat(-G/(H + obj.tree_regularization), size(y, 1), 1);
      predictor = obj.tree_predictor();
      if isprop(predictor, 'weak_models')
        predictor = predictor.setUseModel(modelID);
      end
      predictor = predictor.trainModel(X, w);
    end
    
    function pruneRecursive(obj, X, y, iNode)
      if obj.tree_nodes(iNode).left > 0 && obj.tree_nodes(iNode).right > 0 ...
        && obj.tree_nodes(obj.tree_nodes(iNode).left) == 0 ...
        && obj.tree_nodes(obj.tree_nodes(iNode).right) == 0
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        
        left = struct('idx', obj.tree_nodes(iNode).splitter(X), 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~left.idx, 'iNode', obj.tree_nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        G = sum(y(:, 1));
        H = sum(y(:, 2));
        w = repmat(-G/(H + obj.tree_regularization), size(X, 1), 1);
        objective = obj.objectiveFunc(w, yPred);
        
        current = struct;
        current.X = [obj.tree_nodes(left.iNode).X; obj.tree_nodes(right.iNode).X];
        current.y = [obj.tree_nodes(left.iNode).y; obj.tree_nodes(right.iNode).y];
        current.xMean = mean(X);
        predictorNew = obj.predictorFunc(current.xMean);
        predictorNew.trainModel(current.X, current.y, current.xMean, 0);
        
        objectiveNew = obj.objectiveFunc(w, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          % remove children
          obj.tree_nodes(left.iNode) = TreeModel.nodeTemplate;
          obj.tree_nodes(right.iNode) = TreeModel.nodeTemplate;
          obj.tree_nodes(iNode).left = 0;
          obj.tree_nodes(iNode).right = 0;
          obj.tree_nodes(iNode).predictor = predictorNew;
          obj.tree_nodes(iNode).X = X;
          obj.tree_nodes(iNode).y = y;
        end
      end
    end
  end
  
end