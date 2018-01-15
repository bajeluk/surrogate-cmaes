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
    tree_nodes         % list of nodes
    tree_nNodes        % nodes count
    tree_minGain       % minimum gain to split
    tree_minLeafSize   % minimum number of data examples in leaf
    tree_minParentSize % minimum number of data examples in parent
    tree_maxDepth      % maximum depth of the tree
    tree_splitFunc     % cell-array of split functions
    tree_splits        % generators for split functions
    tree_splitGainFunc % function defining splitGain
    tree_splitGain     % evaluator for split functions
    tree_predictorFunc % function which creates a model in leaf
    tree_predictor     % predictor with appropriate settings
    tree_growFull      % grows a full tree then prunes, otherwise prunes during splitting
    tree_lossFunc      % loss function used for pruning
    tree_fuzziness     % use fuzzy splits (range [0,1])
  end
  
  methods
    function obj = TreeModel(modelOptions)
      % constructor
      obj = obj@WeakModel(modelOptions);
      obj.tree_minGain = defopts(modelOptions, 'tree_minGain', 1e-2);
      obj.tree_minLeafSize = defopts(modelOptions, 'tree_minLeafSize', 2);
      obj.tree_minParentSize = defopts(modelOptions, 'tree_minParentSize', 10);
      obj.tree_maxDepth = defopts(modelOptions, 'tree_maxDepth', inf);
      obj.tree_growFull = defopts(modelOptions, 'tree_growFull', false);
      obj.tree_lossFunc = defopts(modelOptions, 'tree_lossFunc', @mseLossFunc);
      obj.tree_predictorFunc = defopts(modelOptions, 'tree_predictorFunc', ...
        @ConstantModel);
      obj.tree_predictor = obj.tree_predictorFunc(modelOptions);
      obj.tree_splitFunc = defopts(modelOptions, 'tree_splitFunc', ...
        {@AxisSplit}); 
      if ~iscell(obj.tree_splitFunc)
        obj.tree_splitFunc = {obj.tree_splitFunc};
      end
      obj.tree_splits = cellfun(@(x) x(modelOptions), obj.tree_splitFunc, 'UniformOutput', false);
      obj.tree_splitGainFunc = defopts(modelOptions, 'tree_splitGainFunc', ...
        @MSESplitGain);
      obj.tree_splitGain = obj.tree_splitGainFunc(modelOptions);
      obj.tree_fuzziness = defopts(modelOptions, 'tree_fuzziness', 0);
    end

    function obj = trainModel(obj, X, y)
      % train the model based on the data (X,y)      
      obj.tree_nNodes = 0;
      initialSize = min(2^obj.tree_maxDepth, 2 * round(size(X,1) / obj.tree_minLeafSize));
      obj.tree_nodes = repmat(TreeModel.nodeTemplate, initialSize, 1);
      iNodeRoot = obj.addNode();
      obj.trainModelRecursive(X, y, iNodeRoot, 0);
    end

    function [yPred, sd2, ci] = modelPredict(obj, X)
      % predicts the function values in new points X
      if obj.tree_fuzziness == 0
        [yPred, sd2] = obj.modelPredictRecursive(X, 1);
      else
        [yPred, sd2] = obj.modelPredictFuzzyRecursive(X, 1);
      end
      if nargout >= 3
        ci = varToConfidence(yPred, sd2);
      end
    end
    
    function N = getMinTrainPoints(obj, dim)
    % returns minimal number of points necessary to train the model
      N = obj.tree_predictor.getMinTrainPoints(dim) + 1;
    end
    
    function obj = prune(obj, X, y)
      obj.pruneRecursive(X, y, 1);
    end
  end
  
  methods (Access = protected)
    function trainModelRecursive(obj, X, y, iNode, depth)
      [n, dim] = size(X);
      % count actual minimal leaf size according to data dimension
      minLeafSize = max(obj.tree_minLeafSize, obj.getMinTrainPoints(dim));
      if depth < obj.tree_maxDepth ...
          && n >= obj.tree_minParentSize ...
          && n >= 2 * minLeafSize ...
          && size(unique(y, 'rows'), 1) >= 2
        best = struct('gain', -inf);
        splitGain = obj.tree_splitGain.reset(X, y);
        for iSplit = 1:size(obj.tree_splits, 2)
          split = obj.tree_splits{iSplit}.reset(X, y);
          candidate = split.get(splitGain);
          idx = candidate.splitter(X) <= 0.5;
          if sum(idx) >= minLeafSize && sum(~idx) >= minLeafSize
            if candidate.gain > best.gain
              best = candidate;
            end
          end
        end
        if best.gain > -inf && (best.gain >= obj.tree_minGain || obj.tree_growFull)
          idx = best.splitter(X) <= 0.5;
          
          left = struct('idx', idx, 'iNode', obj.addNode());
          left.X = X(left.idx, :);
          left.y = y(left.idx, :);
          obj.trainModelRecursive(left.X, left.y, left.iNode, depth+1);
          
          right = struct('idx', ~idx, 'iNode', obj.addNode());
          right.X = X(right.idx, :);
          right.y = y(right.idx, :);
          obj.trainModelRecursive(right.X, right.y, right.iNode, depth+1);
          
          if best.gain < obj.tree_minGain && obj.tree_nodes(left.iNode).leaf && obj.tree_nodes(right.iNode).leaf
            % prune
            obj.tree_nodes(left.iNode) = obj.nodeTemplate;
            obj.tree_nodes(right.iNode) = obj.nodeTemplate;
          else
            % keep it
            obj.tree_nodes(iNode).leaf = false;
            obj.tree_nodes(iNode).splitter = best.splitter;
            obj.tree_nodes(iNode).left = left.iNode;
            obj.tree_nodes(iNode).right = right.iNode;
            
            obj.tree_nodes(left.iNode).parent = iNode;
            if obj.tree_nodes(left.iNode).leaf
              obj.tree_nodes(left.iNode).predictor = obj.trainPredictor(left.X, left.y);
              %obj.tree_nodes(left.iNode).X = left.X;
              %obj.tree_nodes(left.iNode).y = left.y;
            end
            
            obj.tree_nodes(right.iNode).parent = iNode;
            if obj.tree_nodes(right.iNode).leaf
              obj.tree_nodes(right.iNode).predictor = obj.trainPredictor(right.X, right.y);
              %obj.tree_nodes(right.iNode).X = right.X;
              %obj.tree_nodes(right.iNode).y = right.y;
            end
          end
        end
      end
      if iNode == 1 && obj.tree_nodes(iNode).leaf
        % predictor in root
        obj.tree_nodes(iNode).predictor = obj.trainPredictor(X, y);
        %obj.tree_nodes(iNode).X = X;
        %obj.tree_nodes(iNode).y = y;
      end
    end
    
    function predictor = trainPredictor(obj, X, y)
      predictor = obj.tree_predictor();
      predictor = predictor.trainModel(X, y);
    end
    
    function [yPred, sd2] = modelPredictRecursive(obj, X, iNode)
      if obj.tree_nodes(iNode).leaf
        [yPred, sd2] = obj.tree_nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        idx = obj.tree_nodes(iNode).splitter(X) <= 0.5;
        
        left = struct('idx', idx, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', ~idx, 'iNode', obj.tree_nodes(iNode).right);
        if any(right.idx)
          [yPred(right.idx), sd2(right.idx)] = obj.modelPredictRecursive(X(right.idx, :), right.iNode);
        end
      end
    end
    
    function [yPred, sd2] = modelPredictFuzzyRecursive(obj, X, iNode)
      if obj.tree_nodes(iNode).leaf
        [yPred, sd2] = obj.tree_nodes(iNode).predictor.modelPredict(X);
      else
        yPred = zeros(size(X, 1), 1);
        sd2 = zeros(size(X, 1), 1);
        p = obj.tree_nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.tree_fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.tree_fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx), sd2(left.idx)] = obj.modelPredictFuzzyRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.tree_nodes(iNode).right);
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
      if obj.tree_nodes(iNode).left > 0 && obj.tree_nodes(iNode).right > 0 ...
        && obj.tree_nodes(obj.tree_nodes(iNode).left) == 0 ...
        && obj.tree_nodes(obj.tree_nodes(iNode).right) == 0
        % internal node having two leaves
        % consider making this node a leaf
        yPred = zeros(size(X, 1), 1);
        p = obj.tree_nodes(iNode).splitter(X);
        pTresholdLeft = 0.5 - 0.5 * obj.tree_fuzziness;
        pTresholdRight = 0.5 + 0.5 * obj.tree_fuzziness;
        
        left = struct('idx', p <= pTresholdLeft, 'iNode', obj.tree_nodes(iNode).left);
        if any(left.idx)
          [yPred(left.idx)] = obj.modelPredictRecursive(X(left.idx, :), left.iNode);
        end
        
        right = struct('idx', p > pTresholdRight, 'iNode', obj.tree_nodes(iNode).right);
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
        
        objective = obj.tree_growFull(y, yPred);
        
        current = struct;
        current.X = [obj.tree_nodes(left.iNode).X; obj.tree_nodes(right.iNode).X];
        current.y = [obj.tree_nodes(left.iNode).y; obj.tree_nodes(right.iNode).y];
        predictorNew = obj.tree_predictor();
        predictorNew = predictorNew.trainModel(current.X, current.y);
        objectiveNew = obj.tree_growFull(y, predictorNew.modelPredict(X));
        
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
      if obj.tree_nodes(iNode).left == 0 && obj.tree_nodes(iNode).right == 0
        % try to replace complex model in leaf with constant model
        objective = obj.tree_growFull(y, obj.tree_nodes(iNode).predictor.modelPredict(X));
        
        predictorNew = ConstantModel(struct);
        predictorNew = predictorNew.trainModel(obj.tree_nodes(iNode).X, obj.tree_nodes(iNode).y);
        objectiveNew = obj.tree_growFull(y, predictorNew.modelPredict(X));
        
        if objectiveNew <= objective
          obj.tree_nodes(iNode).predictor = predictorNew;
        end
      end
    end
    
    function iNode = addNode(obj)
      obj.tree_nNodes = obj.tree_nNodes + 1;
      if size(obj.tree_nodes, 1) < obj.tree_nNodes
        obj.children(2 * obj.tree_nNodes, 1) = TreeModel.nodeTemplate;
      end
      iNode = obj.tree_nNodes;
    end
  end
  
end