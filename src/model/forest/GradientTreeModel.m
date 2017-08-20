classdef GradientTreeModel < TreeModel  
  properties %(Access = protected)
    regularization % Ridge regularization parameter
  end
  
  methods
    function obj = GradientTreeModel(modelOptions, xMean)
      % constructor
      regularization = defopts(modelOptions, 'regularization', 0);
      modelOptions.splitGain = GradientSplitGain(regularization);
      modelOptions.predictorFunc = @(xMean) ConstantModel(struct, xMean);
      modelOptions.objectiveFunc = @immse;
      obj = obj@TreeModel(modelOptions, xMean);
      
      % specific model options
      obj.regularization = defopts(modelOptions, 'regularization', 0);
    end
  end
  
  methods (Access = private)
    function trainModelRecursive(obj, X, y, iNode, depth)
      [N, D] = size(X);
      if depth < obj.maxDepth && N >= obj.minParentSize && N >= 2 * obj.minLeafSize && length(unique(y)) >= 2
        best = struct('gain', -inf);
        obj.splitGain.reset(X, y);
        for iSplit = 1:size(obj.splits, 2)
          split = obj.splits{iSplit};
          split.reset(X, y);
          candidate = split.get();
          idx = candidate.splitter(X);
          if sum(idx) >= obj.minLeafSize && sum(~idx) >= obj.minLeafSize
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
        if best.gain >= obj.minGain
          obj.nodes(iNode).splitter = best.splitter;
          idx = best.splitter(X);
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          obj.nodes(iNode).left = left.iNode;
          obj.nodes(left.iNode).parent = iNode;
          obj.trainModelRecursive(X(left.idx, :), y(left.idx, :), left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          obj.nodes(iNode).right = right.iNode;
          obj.nodes(right.iNode).parent = iNode;
          obj.trainModelRecursive(X(right.idx, :), y(right.idx, :), right.iNode, depth+1);
          
          return;
        end
      end
      if isempty(obj.nodes(iNode).predictor)
        G = sum(obj.y(:, 1));
        H = sum(obj.y(:, 2));
        w = repmat(-G/(H + obj.regularization), size(X, 1), 1);
        obj.nodes(iNode).predictor = obj.predictorFunc(mean(X));
        obj.nodes(iNode).predictor.trainModel(X, w, mean(X), 0);
        obj.nodes(iNode).X = X;
        obj.nodes(iNode).y = y;
      end
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