classdef GradientTreeModel < TreeModel  
  properties %(Access = protected)
    regularization % Ridge regularization parameter
  end
  
  methods
    function obj = GradientTreeModel(modelOptions)
      % constructor
      splitGainOptions.weightedGain = false;
      splitGainOptions.regularization = defopts(modelOptions, 'regularization', 0);
      modelOptions.splitGain = GradientSplitGain(splitGainOptions);
      modelOptions.predictorFunc = @() ConstantModel(struct);
      modelOptions.lossFunc = @immse;
      obj = obj@TreeModel(modelOptions);
      
      % specific model options
      obj.regularization = defopts(modelOptions, 'regularization', 0);
    end
  end
  
  methods (Access = protected)
    function predictor = trainPredictor(obj, X, y)
      G = sum(y(:, 1));
      H = sum(y(:, 2));
      w = repmat(-G/(H + obj.regularization), size(y, 1), 1);
      predictor = obj.predictorFunc();
      predictor = predictor.trainModel(X, w);
    end
    
    function pruneRecursive(obj, X, y, iNode)
      if obj.nodes(iNode).left > 0 && obj.nodes(iNode).right > 0 ...
        && obj.nodes(obj.nodes(iNode).left) == 0 ...
        && obj.nodes(obj.nodes(iNode).right) == 0
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        
        left = struct('idx', obj.nodes(iNode).splitter(X), 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~left.idx, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        G = sum(y(:, 1));
        H = sum(y(:, 2));
        w = repmat(-G/(H + obj.regularization), size(X, 1), 1);
        objective = obj.objectiveFunc(w, yPred);
        
        current = struct;
        current.X = [obj.nodes(left.iNode).X; obj.nodes(right.iNode).X];
        current.y = [obj.nodes(left.iNode).y; obj.nodes(right.iNode).y];
        current.xMean = mean(X);
        predictorNew = obj.predictorFunc(current.xMean);
        predictorNew.trainModel(current.X, current.y, current.xMean, 0);
        
        objectiveNew = obj.objectiveFunc(w, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          % remove children
          obj.nodes(left.iNode) = TreeModel.nodeTemplate;
          obj.nodes(right.iNode) = TreeModel.nodeTemplate;
          obj.nodes(iNode).left = 0;
          obj.nodes(iNode).right = 0;
          obj.nodes(iNode).predictor = predictorNew;
          obj.nodes(iNode).X = X;
          obj.nodes(iNode).y = y;
        end
      end
    end
  end
  
end