classdef TreeModel < WeakModel
  properties (Constant, Access = protected)
    nodeTemplate = struct(... % template for nodes
        'leaf', true, ...
        'parent', 0, ...
        'left', 0, ...
        'right', 0, ...
        'splitter', [], ...
        'predictor', [], ...
        'X', [], ...
        'y', []);
  end
  
  properties %(Access = protected)
    nodes % list of nodes
    nNodes % nodes count
    minGain % minimum gain to split
    minLeafSize % minimum number of data examples in leaf
    minParentSize % minimum number of data examples in parent
    maxDepth % maximum depth of the tree
    splits % generators for split functions
    splitGain % evaluator for split functions
    predictorFunc % function which creates a model in leaf
    growFull % grows a full tree then prunes, otherwise prunes during splitting
    lossFunc % loss function used for pruning
    fuzziness % use fuzzy splits (range [0,1])
  end
  
  methods
    function obj = TreeModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      obj.minGain = defopts(modelOptions, 'minGain', 1e-2);
      obj.minLeafSize = defopts(modelOptions, 'minLeafSize', 2);
      obj.minParentSize = defopts(modelOptions, 'minParentSize', 10);
      obj.maxDepth = defopts(modelOptions, 'maxDepth', inf);
      obj.growFull = defopts(modelOptions, 'growFull', false);
      obj.lossFunc = defopts(modelOptions, 'lossFunc', @immse);
      obj.predictorFunc = defopts(modelOptions, 'predictorFunc', ...
        @() ConstantModel(struct));
      obj.splits = defopts(modelOptions, 'splits', ...
        {AxisSplit(struct)});
      obj.splitGain = defopts(modelOptions, 'splitGain', ...
        MSESplitGain(struct));
      obj.fuzziness = defopts(modelOptions, 'fuzziness', 0);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)      
      obj.nNodes = 0;
      initialSize = min(2^obj.maxDepth, 2 * round(size(X,1) / obj.minLeafSize));
      obj.nodes = repmat(TreeModel.nodeTemplate, initialSize, 1);
      iNodeRoot = obj.addNode();
      obj.trainModelRecursive(X, y, iNodeRoot, 0);
    end

    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      if obj.fuzziness == 0
        [yPred, sd2] = obj.modelPredictRecursive(X, 1);
      else
        [yPred, sd2] = obj.modelPredictFuzzyRecursive(X, 1);
      end
      if nargout >= 3
        ci = varToConfidence(yPred, sd2);
      end
    end
    
    function obj = prune(obj, X, y)
      obj.pruneRecursive(X, y, 1);
    end
  end
  
  methods (Access = protected)
    function trainModelRecursive(obj, X, y, iNode, depth)
      n = size(X, 1);
      if depth < obj.maxDepth ...
          && n >= obj.minParentSize && n >= 2 * obj.minLeafSize ...
          && size(unique(y, 'rows'), 1) >= 2
        best = struct('gain', -inf);
        splitGain = obj.splitGain.reset(X, y);
        for iSplit = 1:size(obj.splits, 2)
          split = obj.splits{iSplit}.reset(X, y);
          candidate = split.get(splitGain);
          idx = candidate.splitter(X) <= 0.5;
          if sum(idx) >= obj.minLeafSize && sum(~idx) >= obj.minLeafSize
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
        if best.gain > -inf && (best.gain >= obj.minGain || obj.growFull)
          idx = best.splitter(X) <= 0.5;
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          left.X = X(left.idx, :);
          left.y = y(left.idx, :);
          obj.trainModelRecursive(left.X, left.y, left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          right.X = X(right.idx, :);
          right.y = y(right.idx, :);
          obj.trainModelRecursive(right.X, right.y, right.iNode, depth+1);
          
          if best.gain < obj.minGain && obj.nodes(left.iNode).leaf && obj.nodes(right.iNode).leaf
            % prune
            obj.nodes(left.iNode) = obj.nodeTemplate;
            obj.nodes(right.iNode) = obj.nodeTemplate;
          else
            % keep it
            obj.nodes(iNode).leaf = false;
            obj.nodes(iNode).splitter = best.splitter;
            obj.nodes(iNode).left = left.iNode;
            obj.nodes(iNode).right = right.iNode;
            
            obj.nodes(left.iNode).parent = iNode;
            if obj.nodes(left.iNode).leaf
              obj.nodes(left.iNode).predictor = obj.trainPredictor(left.X, left.y);
              %obj.nodes(left.iNode).X = left.X;
              %obj.nodes(left.iNode).y = left.y;
            end
            
            obj.nodes(right.iNode).parent = iNode;
            if obj.nodes(right.iNode).leaf
              obj.nodes(right.iNode).predictor = obj.trainPredictor(right.X, right.y);
              %obj.nodes(right.iNode).X = right.X;
              %obj.nodes(right.iNode).y = right.y;
            end
          end
        end
      end
      if iNode == 1 && obj.nodes(iNode).leaf
        % predictor in root
        obj.nodes(iNode).predictor = obj.trainPredictor(X, y);
        %obj.nodes(iNode).X = X;
        %obj.nodes(iNode).y = y;
      end
    end
    
    function predictor = trainPredictor(obj, X, y)
      predictor = obj.predictorFunc();
      predictor = predictor.trainModel(X, y);
    end
    
    function [yPred, sd2] = modelPredictRecursive(obj, X, iNode)
      if obj.nodes(iNode).leaf
        [yPred, sd2] = obj.nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        idx = obj.nodes(iNode).splitter(X) <= 0.5;
        
        left = struct('idx', idx, 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~idx, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
      end
    end
    
    function [yPred, sd2] = modelPredictFuzzyRecursive(obj, X, iNode)
      if obj.nodes(iNode).leaf
        [yPred, sd2] = obj.nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        p = obj.nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictFuzzyRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictFuzzyRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > pTresholdLeft, p <= pTresholdRight));
        if any(both.idx)
          both.X = X(both.idx, :);
          pRight = p(both.idx);
          pLeft = 1 - pRight;
          [yLeft, sd2Left] = obj.modelPredictFuzzyRecursive(both.X, left.iNode);
          [yRight, sd2Right] = obj.modelPredictFuzzyRecursive(both.X, right.iNode);
          [yPred(both.idx)] = pLeft .* yLeft + pRight .* yRight;
          % TODO can we do this?
          [sd2(both.idx)] = pLeft .* sd2Left + pRight .* sd2Right;
          % [sd2(both.idx)] = pLeft.^2 .* sd2Left + pRight.^2 .* sd2Right;
        end
      end
    end
    
    function pruneRecursive(obj, X, y, iNode)
      if obj.nodes(iNode).left > 0 && obj.nodes(iNode).right > 0 ...
        && obj.nodes(obj.nodes(iNode).left) == 0 ...
        && obj.nodes(obj.nodes(iNode).right) == 0
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        p = obj.nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > pTresholdLeft, p <= pTresholdRight));
        if any(both.idx)
          both.X = X(both.idx, :);
          pRight = p(both.idx);
          pLeft = 1 - pRight;
          [yLeft] = obj.modelPredictFuzzyRecursive(both.X, left.iNode);
          [yRight] = obj.modelPredictFuzzyRecursive(both.X, right.iNode);
          [yPred(both.idx)] = pLeft .* yLeft + pRight .* yRight;
        end
        
        objective = obj.lossFunc(y, yPred);
        
        current = struct;
        current.X = [obj.nodes(left.iNode).X; obj.nodes(right.iNode).X];
        current.y = [obj.nodes(left.iNode).y; obj.nodes(right.iNode).y];
        predictorNew = obj.predictorFunc();
        predictorNew = predictorNew.trainModel(current.X, current.y);
        objectiveNew = obj.lossFunc(y, predictorNew.modelPredict(X));
        
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
      if obj.nodes(iNode).left == 0 && obj.nodes(iNode).right == 0
        % try to replace complex model in leaf with constant model
        objective = obj.lossFunc(y, obj.nodes(iNode).predictor.modelPredict(X));
        
        predictorNew = ConstantModel(struct);
        predictorNew = predictorNew.trainModel(obj.nodes(iNode).X, obj.nodes(iNode).y);
        objectiveNew = obj.lossFunc(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          obj.nodes(iNode).predictor = predictorNew;
        end
      end
    end
    
    function iNode = addNode(obj)
      obj.nNodes = obj.nNodes + 1;
      if size(obj.nodes, 1) < obj.nNodes
        obj.children(2 * obj.nNodes, 1) = TreeModel.nodeTemplate;
      end
      iNode = obj.nNodes;
    end
  end
  
end