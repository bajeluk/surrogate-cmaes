classdef TreeModel < WeakModel
  properties (Constant, Access = private)
    nodeTemplate = struct(... % template for nodes
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
    pruning % grows a full tree then prunes, otherwise prunes during splitting
    lossFunc % loss function used for pruning
    fuzziness % use fuzzy splits
  end
  
  methods
    function obj = TreeModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      obj.minGain = defopts(modelOptions, 'minGain', 1e-2);
      obj.minLeafSize = defopts(modelOptions, 'minLeafSize', 5);
      obj.minParentSize = defopts(modelOptions, 'minParentSize', 10);
      obj.maxDepth = defopts(modelOptions, 'maxDepth', inf);
      obj.pruning = defopts(modelOptions, 'pruning', false);
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
  
  methods (Access = private)
    function trainModelRecursive(obj, X, y, iNode, depth)
      n = size(X, 1);
      if depth < obj.maxDepth && n >= obj.minParentSize && n >= 2 * obj.minLeafSize && length(unique(y)) >= 2
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
        if best.gain >= obj.minGain
          obj.nodes(iNode).splitter = best.splitter;
          idx = best.splitter(X) <= 0.5;
          
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
        obj.nodes(iNode).predictor = obj.predictorFunc();
        obj.nodes(iNode).predictor = obj.nodes(iNode).predictor.trainModel(X, y);
        obj.nodes(iNode).X = X;
        obj.nodes(iNode).y = y;
      end
    end
    
    function [yPred, sd2] = modelPredictRecursive(obj, X, iNode)
      if isempty(obj.nodes(iNode).splitter)
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
      if isempty(obj.nodes(iNode).splitter)
        [yPred, sd2] = obj.nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        p = obj.nodes(iNode).splitter(X);
        
        left = struct('idx', p <= 0.5 - obj.fuzziness, 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictFuzzyRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > 0.5 + obj.fuzziness, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictFuzzyRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > 0.5 - obj.fuzziness, p <= 0.5 + obj.fuzziness));
        if any(both.idx)
          XBoth = X(both.idx, :);
          pRight = p(both.idx);
          [yLeft, sd2Left] = obj.modelPredictFuzzyRecursive(XBoth, left.iNode);
          [yRight, sd2Right] = obj.modelPredictFuzzyRecursive(XBoth, right.iNode);
          [yPred(both.idx)] = (1-pRight) .* yLeft + pRight .* yRight;
          % TODO can we do this?
          [sd2(both.idx)] = (1-pRight) .* sd2Left + pRight .* sd2Right;
          % [sd2(both.idx)] = (1-p(both.idx)).^2 .* sd2Left + p(both.idx).^2 .* sd2Right;
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
        
        left = struct('idx', p <= 0.5 - obj.fuzziness, 'iNode', obj.nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > 0.5 + obj.fuzziness, 'iNode', obj.nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
        
        both = struct('idx', and(p > 0.5 - obj.fuzziness, p <= 0.5 + obj.fuzziness));
        if any(both.idx)
          XBoth = X(both.idx, :);
          pRight = p(both.idx);
          [yLeft] = obj.modelPredictFuzzyRecursive(XBoth, left.iNode);
          [yRight] = obj.modelPredictFuzzyRecursive(XBoth, right.iNode);
          [yPred(both.idx)] = (1-pRight) .* yLeft + pRight .* yRight;
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